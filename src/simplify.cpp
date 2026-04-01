#include "simplify.h"

#include <cmath>
#include <optional>
#include <queue>
#include <vector>

namespace {

void append_vertex_if_distinct(std::vector<Point>& vertices, const Point& p) {
    if (vertices.empty() || !point_equal(vertices.back(), p)) {
        vertices.push_back(p);
    }
}

void normalize_ring_vertices(std::vector<Point>& vertices) {
    if (vertices.empty()) {
        return;
    }

    std::vector<Point> cleaned;
    cleaned.reserve(vertices.size());
    for (const Point& p : vertices) {
        append_vertex_if_distinct(cleaned, p);
    }

    if (cleaned.size() >= 2 && point_equal(cleaned.front(), cleaned.back())) {
        cleaned.pop_back();
    }

    vertices.swap(cleaned);
}

Ring collapsed_ring(const Ring& ring, int iA, int iB, int iC, int iD, const Point& E) {
    (void)iA;
    (void)iD;

    Ring out;
    out.ring_id = ring.ring_id;

    const int n = static_cast<int>(ring.vertices.size());
    out.vertices.reserve(n - 1);

    for (int i = 0; i < n; ++i) {
        if (i == iB) {
            append_vertex_if_distinct(out.vertices, E);
        }
        if (i == iB || i == iC) {
            continue;
        }
        append_vertex_if_distinct(out.vertices, ring.vertices[i]);
    }

    normalize_ring_vertices(out.vertices);
    return out;
}

bool candidate_preserves_minimum_ring_size(const Ring& ring) {
    return static_cast<int>(ring.vertices.size()) >= 5;
}

bool same_directed_edge(int i, int j, int u, int v) {
    return i == u && j == v;
}

bool segment_intersects_ring_edges_except_window(
    const Segment& s,
    const Ring& ring,
    int iA, int iB, int iC, int iD,
    bool is_AE)
{
    const int n = static_cast<int>(ring.vertices.size());
    const int iPrevA = (iA - 1 + n) % n;
    const int iNextD = (iD + 1) % n;

    for (int i = 0; i < n; ++i) {
        const int j = (i + 1) % n;

        if (same_directed_edge(i, j, iA, iB) ||
            same_directed_edge(i, j, iB, iC) ||
            same_directed_edge(i, j, iC, iD)) {
            continue;
        }

        if (is_AE && same_directed_edge(i, j, iPrevA, iA)) {
            continue;
        }
        if (!is_AE && same_directed_edge(i, j, iD, iNextD)) {
            continue;
        }

        const Point& p = ring.vertices[i];
        const Point& q = ring.vertices[j];

        if (segments_intersect(s.a, s.b, p, q)) {
            return true;
        }
    }
    return false;
}

bool segment_intersects_ring_edges(const Segment& s, const Ring& ring) {
    const int n = static_cast<int>(ring.vertices.size());
    for (int i = 0; i < n; ++i) {
        const Point& p = ring.vertices[i];
        const Point& q = ring.vertices[(i + 1) % n];

        if (segments_intersect(s.a, s.b, p, q)) {
            return true;
        }
    }
    return false;
}

bool candidate_topology_valid(
    const Polygon& poly,
    int ring_index,
    int iA, int iB, int iC, int iD,
    const Point& E)
{
    const Ring& ring = poly.rings[ring_index];
    const Point& A = ring.vertices[iA];
    const Point& D = ring.vertices[iD];

    Segment AE{A, E};
    Segment ED{E, D};

    if (!point_equal(A, E)) {
        if (segment_intersects_ring_edges_except_window(AE, ring, iA, iB, iC, iD, true)) return false;
        for (int r = 0; r < static_cast<int>(poly.rings.size()); ++r) {
            if (r == ring_index) continue;
            if (segment_intersects_ring_edges(AE, poly.rings[r])) return false;
        }
    }

    if (!point_equal(E, D)) {
        if (segment_intersects_ring_edges_except_window(ED, ring, iA, iB, iC, iD, false)) return false;
        for (int r = 0; r < static_cast<int>(poly.rings.size()); ++r) {
            if (r == ring_index) continue;
            if (segment_intersects_ring_edges(ED, poly.rings[r])) return false;
        }
    }

    Polygon test = poly;
    test.rings[ring_index] = collapsed_ring(ring, iA, iB, iC, iD, E);
    if (static_cast<int>(test.rings[ring_index].vertices.size()) < 3) {
        return false;
    }
    return polygon_topology_valid(test);
}

std::optional<Candidate> make_candidate(
    const Polygon& poly,
    int ring_index,
    int iA)
{
    const Ring& ring = poly.rings[ring_index];
    const int n = static_cast<int>(ring.vertices.size());

    if (!candidate_preserves_minimum_ring_size(ring)) return std::nullopt;
    if (n < 4) return std::nullopt;

    const int iB = (iA + 1) % n;
    const int iC = (iA + 2) % n;
    const int iD = (iA + 3) % n;

    const Point& A = ring.vertices[iA];
    const Point& B = ring.vertices[iB];
    const Point& C = ring.vertices[iC];
    const Point& D = ring.vertices[iD];

    ApscPlacement apsc = apsc_placement_and_displacement(A, B, C, D);
    if (!apsc.valid) return std::nullopt;

    Ring new_ring = collapsed_ring(ring, iA, iB, iC, iD, apsc.E);
    if (static_cast<int>(new_ring.vertices.size()) < 3) return std::nullopt;
    if (!candidate_topology_valid(poly, ring_index, iA, iB, iC, iD, apsc.E)) return std::nullopt;

    Candidate c;
    c.ring_index = ring_index;
    c.iA = iA;
    c.iB = iB;
    c.iC = iC;
    c.iD = iD;
    c.E = apsc.E;
    c.displacement = apsc.displacement;
    return c;
}

struct PQItem {
    Candidate cand;
    std::size_t generation;
};

struct PQCompare {
    bool operator()(const PQItem& lhs, const PQItem& rhs) const {
        if (std::fabs(lhs.cand.displacement - rhs.cand.displacement) > EPS) {
            return lhs.cand.displacement > rhs.cand.displacement;
        }
        if (lhs.cand.ring_index != rhs.cand.ring_index) {
            return lhs.cand.ring_index > rhs.cand.ring_index;
        }
        return lhs.cand.iA > rhs.cand.iA;
    }
};

void apply_candidate(Polygon& poly, const Candidate& cand) {
    Ring& ring = poly.rings[cand.ring_index];
    ring = collapsed_ring(ring, cand.iA, cand.iB, cand.iC, cand.iD, cand.E);
}

} // namespace

Polygon simplify_polygon(const Polygon& input, int target_vertices, double& total_areal_displacement) {
    Polygon poly = input;
    total_areal_displacement = 0.0;

    if (target_vertices < 3) {
        target_vertices = 3;
    }

    std::vector<std::size_t> ring_generation(poly.rings.size(), 0);
    std::priority_queue<PQItem, std::vector<PQItem>, PQCompare> pq;

    auto push_ring_candidates = [&](int ring_index) {
        const Ring& ring = poly.rings[ring_index];
        const int n = static_cast<int>(ring.vertices.size());
        for (int i = 0; i < n; ++i) {
            auto cand = make_candidate(poly, ring_index, i);
            if (cand.has_value()) {
                pq.push(PQItem{*cand, ring_generation[ring_index]});
            }
        }
    };

    for (int r = 0; r < static_cast<int>(poly.rings.size()); ++r) {
        push_ring_candidates(r);
    }

    while (total_vertices(poly) > target_vertices) {
        bool applied = false;

        while (!pq.empty()) {
            const PQItem top = pq.top();
            pq.pop();

            if (top.generation != ring_generation[top.cand.ring_index]) {
                continue;
            }

            auto fresh = make_candidate(poly, top.cand.ring_index, top.cand.iA);
            if (!fresh.has_value()) {
                continue;
            }

            apply_candidate(poly, *fresh);
            total_areal_displacement += fresh->displacement;
            ++ring_generation[top.cand.ring_index];
            push_ring_candidates(top.cand.ring_index);

            applied = true;
            break;
        }

        if (!applied) {
            break;
        }
    }

    return poly;
}
