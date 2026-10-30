#include "common/zhangProductPublicationBoundary.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

static int assertions = 0;
static void check(bool valid, const char* message) {
    ++assertions;
    if (!valid) throw std::runtime_error(message);
}

// Executes the SAME boundary used by ppp_ambres.cpp and actual tiny directory
// renames. State operands are deliberately an independent stand-in, not claims
// of exercising PEA's real scientific/runtime maps or full-Q writer.
int main(int argc, char** argv) {
    try {
        if (argc > 2) throw std::runtime_error("At most one fresh test output directory is accepted");
        const std::filesystem::path base = argc == 2 ? std::filesystem::path(argv[1]) :
            std::filesystem::temp_directory_path() / ("ginan_publication_boundary_" +
                std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        if (std::filesystem::exists(base)) throw std::runtime_error("Refuse reusing fault-injection directory");
        std::filesystem::create_directories(base);
        for (int stage = 0; stage < 5; ++stage) {
            const auto dir = base / ("case_" + std::to_string(stage));
            const auto pending = dir / "epoch.pending";
            const auto committed = dir / "epoch";
            std::filesystem::create_directories(pending);
            std::ofstream(pending / "COMMIT.json") << "{\"fixture_only\":true}\n";
            std::map<std::string, int> history{{"gamma", 31}, {"counter", 0}};
            const auto before = history;
            int rollbacks = 0;
            bool caught = false;
            try {
                auto rollback = [&]() { ++rollbacks; history = before; };
                ZhangProductPublicationBoundary boundary(rollback);
                check(!boundary.published(), "Boundary began published");
                history = {{"gamma", 0}, {"counter", 1}};
                if (stage == 0) throw std::runtime_error("certificate/full-Q preflight failure");
                if (stage == 1) {
                    std::filesystem::create_directory(committed);
                    std::ofstream(committed / "occupied") << "not a committed fixture\n";
                }
                // Stage 4 abandons the private scope without publishing.
                if (stage != 4) {
                    std::filesystem::rename(pending, committed);
                    boundary.markPublished(); // noexcept; immediate after rename.
                    check(boundary.published(), "Successful rename not marked");
                    if (stage == 2) throw std::runtime_error("post-publication fatal diagnostic failure");
                }
            } catch (const std::exception&) { caught = true; }
            if (stage == 0 || stage == 1 || stage == 4) {
                check(history == before && rollbacks == 1, "Prepublication failure advanced runtime history");
                check(std::filesystem::exists(pending), "Prepublication failure published pending epoch");
                if (stage != 4) check(caught, "Prepublication error became success");
            } else {
                check(history.at("gamma") == 0 && history.at("counter") == 1 && rollbacks == 0,
                      "Postpublication failure restored memory behind visible bundle");
                check(std::filesystem::exists(committed / "COMMIT.json") && !std::filesystem::exists(pending),
                      "Published fixture directory missing");
                check(caught == (stage == 2), "Postpublication error incorrectly reported success");
            }
        }
        std::ofstream out(base / "RESULT.json");
        out << "{\"schema\":\"PUBLICATION_BOUNDARY_REAL_RENAME_FAULT_FIXTURE_V1\",\"passed\":true,\"cases\":5,\"assertions\":"
            << assertions << ",\"production_boundary_helper\":true,\"real_small_directory_rename\":true,\"pea_writer_fault_injection\":false,\"actual_scientific_runtime_maps\":false}\n";
        if (!out) throw std::runtime_error("Failed saving result");
        std::cout << "PASS production boundary helper + real tiny rename faults: " << assertions << " assertions\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
