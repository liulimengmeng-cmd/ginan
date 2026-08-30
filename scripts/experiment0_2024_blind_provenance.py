#!/usr/bin/env python3
"""Leakage-resistant provenance primitives for Experiment 0 2024.

This module deliberately knows nothing about STEC or held-out output formats.
It only records immutable input identities, checks that they did not change
during an evaluation, records the code/runtime identity, and publishes JSON
without replacing an existing result.  The spatial evaluator can therefore
import these functions without giving the provenance layer a reason to open or
interpret any scientific output.

The no-replace writer publishes a fully written temporary file with a hard
link.  Creating that link is atomic and fails if the destination already
exists; unlike ``os.replace``, it cannot overwrite a prior result.  Its trust
boundary is a local filesystem whose parent directories are trusted and are
not concurrently renamed or replaced.  Temporary and final paths are created
in the same directory, so hard-link publication is same-filesystem; filesystems
without atomic hard-link creation are rejected.  Directory fsync is best-effort,
so the guarantee is atomic visibility/no replacement, not power-loss durability.
"""

from __future__ import annotations

import hashlib
import errno
import json
import os
import platform
import shutil
import stat
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np


PROVENANCE_SCHEMA = "GINAN_EXPERIMENT0_2024_BLIND_PROVENANCE_V1"
FILE_RECORD_FIELDS = (
    "path",
    "resolved_path",
    "size",
    "mtime_ns",
    "dev",
    "inode",
    "sha256",
)
_HASH_CHUNK_SIZE = 1024 * 1024

PathLike = str | os.PathLike[str]
FileRecord = dict[str, object]
InputSnapshot = dict[str, FileRecord]


class ProvenanceError(RuntimeError):
    """Base class for provenance contract violations."""


class UnsafePathError(ProvenanceError):
    """Raised for missing, non-regular, or symlinked inputs."""


class InputChangedError(ProvenanceError):
    """Raised when a post-evaluation input snapshot differs from its freeze."""


def _absolute_path(path: PathLike) -> Path:
    return Path(os.path.abspath(os.fspath(path)))


def _is_link_like(path: Path) -> bool:
    if path.is_symlink():
        return True
    try:
        file_attributes = path.lstat().st_file_attributes
    except AttributeError:
        file_attributes = 0
    reparse_point = getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0x400)
    if file_attributes & reparse_point:
        return True
    # Python 3.12 exposes junctions separately on Windows.  Treat them like
    # symlinks when the API is available so a resolved record cannot silently
    # leave the declared directory tree.
    is_junction = getattr(path, "is_junction", None)
    return bool(is_junction and is_junction())


def _reject_link_components(path: Path) -> None:
    """Reject symlinks (and Windows junctions) in every existing component."""

    absolute = _absolute_path(path)
    parts = absolute.parts
    if not parts:
        raise UnsafePathError(f"empty path: {path!s}")
    current = Path(parts[0])
    for part in parts[1:]:
        current /= part
        if os.path.lexists(current) and _is_link_like(current):
            raise UnsafePathError(f"symlink or junction is not allowed: {current}")


def _stat_identity(file_stat: os.stat_result) -> tuple[int, int, int, int]:
    return (
        int(file_stat.st_dev),
        int(file_stat.st_ino),
        int(file_stat.st_size),
        int(file_stat.st_mtime_ns),
    )


def record_regular_file(path: PathLike) -> FileRecord:
    """Return a stable, JSON-ready identity for one regular input file.

    The file is opened without following the final symlink where the platform
    supports ``O_NOFOLLOW``.  Descriptor metadata is checked before and after
    hashing and is also matched to the pathname after hashing.  A file changed
    during capture is rejected instead of producing a mixed record.
    """

    absolute = _absolute_path(path)
    _reject_link_components(absolute)
    try:
        resolved = absolute.resolve(strict=True)
    except FileNotFoundError as exc:
        raise UnsafePathError(f"input does not exist: {absolute}") from exc

    flags = os.O_RDONLY | getattr(os, "O_BINARY", 0)
    flags |= getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(absolute, flags)
    except FileNotFoundError as exc:
        raise UnsafePathError(f"input does not exist: {absolute}") from exc
    except OSError as exc:
        raise UnsafePathError(f"cannot safely open input {absolute}: {exc}") from exc

    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            raise UnsafePathError(f"input is not a regular file: {absolute}")
        digest = hashlib.sha256()
        while True:
            chunk = os.read(descriptor, _HASH_CHUNK_SIZE)
            if not chunk:
                break
            digest.update(chunk)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)

    if _stat_identity(before) != _stat_identity(after):
        raise InputChangedError(f"input changed while hashing: {absolute}")

    # Catch pathname replacement while the descriptor was being read.
    _reject_link_components(absolute)
    try:
        path_after = absolute.stat(follow_symlinks=False)
    except FileNotFoundError as exc:
        raise InputChangedError(f"input disappeared while hashing: {absolute}") from exc
    if not stat.S_ISREG(path_after.st_mode):
        raise UnsafePathError(f"input is no longer a regular file: {absolute}")
    if _stat_identity(after) != _stat_identity(path_after):
        raise InputChangedError(f"input pathname changed while hashing: {absolute}")

    return {
        "path": str(absolute),
        "resolved_path": str(resolved),
        "size": int(after.st_size),
        "mtime_ns": int(after.st_mtime_ns),
        "dev": int(after.st_dev),
        "inode": int(after.st_ino),
        "sha256": digest.hexdigest(),
    }


def capture_input_snapshot(inputs: Mapping[str, PathLike]) -> InputSnapshot:
    """Hash named inputs in deterministic role order."""

    snapshot: InputSnapshot = {}
    for role in sorted(inputs):
        if not isinstance(role, str) or not role:
            raise ProvenanceError("every input role must be a non-empty string")
        snapshot[role] = record_regular_file(inputs[role])
    if not snapshot:
        raise ProvenanceError("at least one input is required")
    return snapshot


def _paths_from_snapshot(
    snapshot: Mapping[str, Mapping[str, object]]
) -> dict[str, str]:
    paths: dict[str, str] = {}
    for role, record in snapshot.items():
        try:
            path = record["path"]
        except KeyError as exc:
            raise ProvenanceError(f"snapshot record {role!r} has no path") from exc
        if not isinstance(path, str) or not path:
            raise ProvenanceError(f"snapshot record {role!r} has an invalid path")
        paths[role] = path
    return paths


def assert_input_snapshot_unchanged(
    before: Mapping[str, Mapping[str, object]],
    inputs: Mapping[str, PathLike] | None = None,
) -> InputSnapshot:
    """Capture a post snapshot and fail if any role or file identity changed.

    ``inputs`` may be omitted when ``before`` was loaded back from JSON; the
    recorded absolute paths are then used.  The unchanged post snapshot is
    returned for inclusion in a final provenance document.
    """

    selected_inputs: Mapping[str, PathLike]
    if inputs is None:
        selected_inputs = _paths_from_snapshot(before)
    else:
        if set(before) != set(inputs):
            raise InputChangedError(
                "input roles changed: "
                f"before={sorted(before)}, after={sorted(inputs)}"
            )
        selected_inputs = inputs

    after = capture_input_snapshot(selected_inputs)
    changes: dict[str, list[str]] = {}
    for role in sorted(before):
        previous = before[role]
        current = after[role]
        changed_fields = [
            field
            for field in FILE_RECORD_FIELDS
            if previous.get(field) != current[field]
        ]
        if changed_fields:
            changes[role] = changed_fields
    if changes:
        raise InputChangedError(
            "input snapshot changed: " + json.dumps(changes, sort_keys=True)
        )
    return after


def _normalise_planned_outputs(
    outputs: Mapping[str, PathLike] | Iterable[PathLike],
) -> dict[str, PathLike]:
    if isinstance(outputs, Mapping):
        named = dict(outputs)
    elif isinstance(outputs, (str, os.PathLike)):
        named = {"output": outputs}
    else:
        named = {f"output_{index}": path for index, path in enumerate(outputs)}
    if not named:
        raise ProvenanceError("at least one planned output is required")
    return named


def assert_planned_outputs_absent(
    outputs: Mapping[str, PathLike] | Iterable[PathLike],
) -> dict[str, str]:
    """Require every planned final path to be absent, including dangling links."""

    planned: dict[str, str] = {}
    seen: dict[str, str] = {}
    for role, path in sorted(_normalise_planned_outputs(outputs).items()):
        if not isinstance(role, str) or not role:
            raise ProvenanceError("every output role must be a non-empty string")
        absolute = _absolute_path(path)
        if os.path.lexists(absolute):
            raise FileExistsError(f"planned output already exists: {absolute}")
        _reject_link_components(absolute.parent)
        resolved = absolute.resolve(strict=False)
        canonical = os.path.normcase(str(resolved))
        if canonical in seen:
            raise ProvenanceError(
                f"planned outputs {seen[canonical]!r} and {role!r} resolve to {resolved}"
            )
        seen[canonical] = role
        planned[role] = str(resolved)
    return planned


def _fsync_directory_best_effort(directory: Path) -> None:
    if os.name != "posix":
        return
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
    try:
        descriptor = os.open(directory, flags)
    except OSError:
        return
    try:
        os.fsync(descriptor)
    except OSError:
        pass
    finally:
        os.close(descriptor)


def _unlink_best_effort(path: Path) -> None:
    try:
        path.unlink()
    except OSError:
        # Once the hard link has been published, failure to remove its private
        # temporary alias must not turn a successful commit into an ambiguous
        # failure.  The dot-prefixed alias is safe to remove later.
        pass


def atomic_write_json_no_replace(path: PathLike, payload: object) -> FileRecord:
    """Atomically publish JSON and fail rather than replace an existing file."""

    destination = _absolute_path(path)
    # Serialize first so an unsupported value cannot leave filesystem debris.
    encoded = (
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")

    if os.path.lexists(destination):
        raise FileExistsError(f"output already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    _reject_link_components(destination.parent)
    if os.path.lexists(destination):
        raise FileExistsError(f"output already exists: {destination}")

    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(encoded)
            stream.flush()
            os.fsync(stream.fileno())
        temporary_record = record_regular_file(temporary)
        published_record = dict(temporary_record)
        published_record["path"] = str(destination)
        published_record["resolved_path"] = str(destination.resolve(strict=False))
        # Hard-link publication has atomic O_EXCL-like semantics on both POSIX
        # and NTFS: EEXIST wins if another process published first.
        try:
            os.link(temporary, destination, follow_symlinks=False)
        except OSError as exc:
            if isinstance(exc, FileExistsError):
                raise
            if exc.errno in {
                errno.EPERM,
                errno.EACCES,
                errno.EXDEV,
                getattr(errno, "ENOTSUP", errno.EINVAL),
                getattr(errno, "EOPNOTSUPP", errno.EINVAL),
            }:
                raise ProvenanceError(
                    "filesystem cannot provide atomic no-replace hard-link "
                    f"publication for {destination}: {exc}"
                ) from exc
            raise
        _fsync_directory_best_effort(destination.parent)
    except BaseException:
        _unlink_best_effort(temporary)
        raise
    else:
        _unlink_best_effort(temporary)
    # No fallible pathname validation follows the os.link commit point.  This
    # identity was measured from the same inode immediately before publication.
    return published_record


def _git(
    executable: str,
    repo: str | os.PathLike[str],
    *arguments: str,
    allow_detached: bool = False,
) -> str | None:
    completed = subprocess.run(
        [executable, "-C", str(repo), *arguments],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
    )
    if allow_detached and completed.returncode == 1:
        return None
    if completed.returncode != 0:
        error = completed.stderr.strip() or completed.stdout.strip()
        raise ProvenanceError(
            f"{executable} {' '.join(arguments)} failed with "
            f"{completed.returncode}: {error}"
        )
    return completed.stdout.strip()


def _wslpath(option: str, path: str | os.PathLike[str]) -> str:
    completed = subprocess.run(
        ["wslpath", option, str(path)],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
    )
    if completed.returncode != 0 or not completed.stdout.strip():
        error = completed.stderr.strip() or completed.stdout.strip()
        raise ProvenanceError(
            f"wslpath {option} failed with {completed.returncode}: {error}"
        )
    return completed.stdout.strip()


def git_identity(repo: PathLike) -> dict[str, object]:
    """Record commit, branch, and cleanliness of tracked files only.

    A worktree created by Windows Git stores a Windows path in its ``.git``
    redirection file.  Native WSL Git interprets that path as a relative Linux
    path and fails.  On WSL only, the function therefore falls back to the
    installed Windows ``git.exe`` and converts paths with ``wslpath``.  The
    recorded repository remains the canonical Linux path, so freeze/evaluate
    identity comparison is stable within the formal WSL runtime.
    """

    requested = _absolute_path(repo)
    executable = "git"
    repository_argument: str = str(requested)
    try:
        top_level_text = _git(
            executable, repository_argument, "rev-parse", "--show-toplevel"
        )
        assert top_level_text is not None
        top_level = Path(top_level_text).resolve(strict=True)
        repository_argument = str(top_level)
    except ProvenanceError as native_error:
        if platform.system() != "Linux" or shutil.which("git.exe") is None:
            raise
        executable = "git.exe"
        repository_argument = _wslpath("-w", requested)
        try:
            top_level_text = _git(
                executable,
                repository_argument,
                "rev-parse",
                "--show-toplevel",
            )
        except ProvenanceError as windows_error:
            raise ProvenanceError(
                f"native Git failed ({native_error}); Windows Git fallback also "
                f"failed ({windows_error})"
            ) from windows_error
        assert top_level_text is not None
        top_level = Path(_wslpath("-u", top_level_text)).resolve(strict=True)
        repository_argument = top_level_text

    if (
        requested.resolve(strict=True) != top_level
        and top_level not in requested.parents
    ):
        raise ProvenanceError(
            f"requested path {requested} is outside Git top level {top_level}"
        )

    head = _git(executable, repository_argument, "rev-parse", "HEAD")
    branch = _git(
        executable,
        repository_argument,
        "symbolic-ref",
        "--quiet",
        "--short",
        "HEAD",
        allow_detached=True,
    )
    status_text = _git(
        executable,
        repository_argument,
        "status",
        "--porcelain=v1",
        "--untracked-files=no",
        "--ignore-submodules=none",
    )
    assert head is not None and status_text is not None
    tracked_status = status_text.splitlines() if status_text else []
    return {
        "repository": str(top_level),
        "backend": executable,
        "head": head,
        "branch": branch,
        "detached": branch is None,
        "tracked_clean": not tracked_status,
        "tracked_status_porcelain_v1": tracked_status,
    }


def runtime_identity() -> dict[str, object]:
    """Return the Python and NumPy runtime identity needed for replay."""

    executable = Path(sys.executable).resolve(strict=True)
    numpy_file = Path(np.__file__).resolve(strict=True) if np.__file__ else None
    return {
        "python": {
            "implementation": platform.python_implementation(),
            "version": platform.python_version(),
            "version_info": list(sys.version_info[:5]),
            "executable": str(executable),
            "compiler": platform.python_compiler(),
        },
        "numpy": {
            "version": np.__version__,
            "module_path": str(numpy_file) if numpy_file else None,
        },
        "platform": platform.platform(),
    }


def begin_provenance(
    repo: PathLike,
    inputs: Mapping[str, PathLike],
    planned_outputs: Mapping[str, PathLike] | Iterable[PathLike],
) -> dict[str, object]:
    """Convenience entry point to call immediately before evaluation."""

    planned = assert_planned_outputs_absent(planned_outputs)
    return {
        "schema": PROVENANCE_SCHEMA,
        "captured_pre_utc": datetime.now(timezone.utc).isoformat(),
        "inputs_pre": capture_input_snapshot(inputs),
        "planned_outputs": planned,
        "git": git_identity(repo),
        "runtime": runtime_identity(),
    }


def complete_provenance(
    provenance: Mapping[str, object],
    inputs: Mapping[str, PathLike] | None = None,
) -> dict[str, object]:
    """Return a completed copy after proving every frozen input is unchanged."""

    if provenance.get("schema") != PROVENANCE_SCHEMA:
        raise ProvenanceError("unsupported or missing provenance schema")
    before = provenance.get("inputs_pre")
    if not isinstance(before, Mapping):
        raise ProvenanceError("provenance document has no inputs_pre mapping")
    completed = dict(provenance)
    completed["inputs_post"] = assert_input_snapshot_unchanged(before, inputs)
    completed["captured_post_utc"] = datetime.now(timezone.utc).isoformat()
    completed["inputs_unchanged"] = True
    return completed
