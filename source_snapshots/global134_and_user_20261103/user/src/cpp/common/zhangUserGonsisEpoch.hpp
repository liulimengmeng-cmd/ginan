#pragma once

struct KFState;
struct KFMeas;

// Bounded, same-epoch, read-only residual evidence for the GONSIS output gate.
// Called only after the authoritative PPP observation transaction commits.
// No observation factor, KF update, NIS or integer history is created here.
void captureZhangUserGonsisFloatObservations(const KFState& state, const KFMeas& measurement);
