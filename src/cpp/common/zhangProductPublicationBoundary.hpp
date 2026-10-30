#pragma once

#include <exception>

/** Runtime rollback boundary around a private product bundle.
 *
 * The callback is referenced, not copied/allocated. A successful filesystem
 * publication MUST be followed immediately by markPublished(), which cannot
 * throw. Before that point an exception or abandoned scope restores the
 * caller's existing runtime snapshots. Afterwards errors remain fatal but
 * must not restore old memory behind an already-visible committed bundle.
 *
 * This is not an output writer, scientific validator or retry mechanism.
 */
template<class Rollback>
class ZhangProductPublicationBoundary
{
    Rollback& rollback_;
    bool published_ = false;

public:
    explicit ZhangProductPublicationBoundary(Rollback& rollback) noexcept
        : rollback_(rollback) {}

    ZhangProductPublicationBoundary(const ZhangProductPublicationBoundary&) = delete;
    ZhangProductPublicationBoundary& operator=(const ZhangProductPublicationBoundary&) = delete;

    void markPublished() noexcept { published_ = true; }
    bool published() const noexcept { return published_; }

    ~ZhangProductPublicationBoundary() noexcept
    {
        if (published_) return;
        try { rollback_(); }
        catch (...)
        {
            // A partial/failed restoration cannot safely continue or replace
            // the active exception with apparent success.
            std::terminate();
        }
    }
};
