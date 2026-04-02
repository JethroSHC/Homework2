#pragma once

#include <optional>
#include "polygon.h"

/**
 * @struct Candidate
 * @brief Represents a possible simplification operation.
 *
 * A candidate replaces vertices B and C in the sequence:
 * A -> B -> C -> D
 * with a new point E, forming:
 * A -> E -> D
*/
struct Candidate {
    int ring_index{-1};
    int iA{-1};
    int iB{-1};
    int iC{-1};
    int iD{-1};
    Point E{};
    double displacement{0.0};

    // stale/valid flags for lazy priority-queue invalidation
    bool valid{true};
    std::size_t generation{0};
};

/**
 * @brief Simplifies a polygon while preserving area and topology.
 *
 * Repeatedly:
 * 1. Generates all valid candidates across all rings
 * 2. Selects the one with smallest displacement
 * 3. Applies the collapse
 *
 * Stops when:
 * - total vertices <= target_vertices
 * - or no valid collapse exists
 *
 * @param input Original polygon
 * @param target_vertices Desired maximum vertex count
 * @param total_areal_displacement Accumulated local areal displacement
 * @return Simplified polygon
 */
Polygon simplify_polygon(const Polygon& input, int target_vertices, double& total_areal_displacement);