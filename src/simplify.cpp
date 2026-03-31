#include "simplify.h"

#include <limits>
#include <optional>

static Ring collapsed_ring(const Ring& ring, int iA, int iB, int iC, int iD, const Point& E) {
    (void)iA;
    (void)iD;

    Ring out;
    out.ring_id = ring.ring_id;

    const int n = static_cast<int>(ring.vertices.size());
    for (int i = 0; i < n; ++i) {
        if (i == iB) {
            out.vertices.push_back(E);
        }
        if (i == iB || i == iC) {
            continue;
        }
        out.vertices.push_back(ring.vertices[i]);
    }

    return out;
}

static bool candidate_preserves_minimum_ring_size(const Ring& ring) {
    return static_cast<int>(ring.vertices.size()) >= 5;
}

static bool candidate_topology_valid(const Polygon& poly, int ring_index, const Ring& new_ring) {
    Polygon test = poly;
    test.rings[ring_index] = new_ring;
    return polygon_topology_valid(test);
}

static std::optional<Candidate> best_candidate_for_polygon(const Polygon& poly) {
    std::optional<Candidate> best;
    double best_disp = std::numeric_limits<double>::infinity();

    for (int r = 0; r < static_cast<int>(poly.rings.size()); ++r) {
        const Ring& ring = poly.rings[r];
        const int n = static_cast<int>(ring.vertices.size());

        if (!candidate_preserves_minimum_ring_size(ring)) continue;

        for (int i = 0; i < n; ++i) {
            const int iA = i;
            const int iB = (i + 1) % n;
            const int iC = (i + 2) % n;
            const int iD = (i + 3) % n;

            const Point& A = ring.vertices[iA];
            const Point& B = ring.vertices[iB];
            const Point& C = ring.vertices[iC];
            const Point& D = ring.vertices[iD];

            ApscPlacement apsc = apsc_placement_and_displacement(A, B, C, D);
            if (!apsc.valid) continue;

            Point E = apsc.E;
            double disp = apsc.displacement;

            Ring new_ring = collapsed_ring(ring, iA, iB, iC, iD, E);

            if (static_cast<int>(new_ring.vertices.size()) < 3) continue;
            if (!candidate_topology_valid(poly, r, new_ring)) continue;

            if (disp < best_disp) {
                best_disp = disp;
                best = Candidate{r, iA, iB, iC, iD, E, disp};
            }
        }
    }

    return best;
}

static void apply_candidate(Polygon& poly, const Candidate& cand) {
    Ring& ring = poly.rings[cand.ring_index];
    ring = collapsed_ring(ring, cand.iA, cand.iB, cand.iC, cand.iD, cand.E);
}

Polygon simplify_polygon(const Polygon& input, int target_vertices, double& total_areal_displacement) {
    Polygon poly = input;
    total_areal_displacement = 0.0;

    if (target_vertices < 3) {
        return poly;
    }

    while (total_vertices(poly) > target_vertices) {
        auto cand = best_candidate_for_polygon(poly);
        if (!cand.has_value()) {
            break;
        }

        apply_candidate(poly, *cand);
        total_areal_displacement += cand->displacement;
    }

    return poly;
}