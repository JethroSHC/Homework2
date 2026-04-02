#include "polygon.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {

void remove_consecutive_duplicate_points(std::vector<Point>& vertices) {
    std::vector<Point> cleaned;
    cleaned.reserve(vertices.size());

    for (const Point& p : vertices) {
        if (cleaned.empty() || !point_equal(cleaned.back(), p)) {
            cleaned.push_back(p);
        }
    }

    if (cleaned.size() >= 2 && point_equal(cleaned.front(), cleaned.back())) {
        cleaned.pop_back();
    }

    vertices.swap(cleaned);
}

bool edges_are_adjacent(int i, int j, int n) {
    if (i == j) return true;
    if ((i + 1) % n == j) return true;
    if ((j + 1) % n == i) return true;
    return false;
}

bool point_on_polygon_boundary(const Point& p, const Polygon& poly) {
    for (const auto& ring : poly.rings) {
        const int n = static_cast<int>(ring.vertices.size());
        for (int i = 0; i < n; ++i) {
            if (on_segment(ring.vertices[i], ring.vertices[(i + 1) % n], p)) {
                return true;
            }
        }
    }
    return false;
}

struct EdgeInfo {
    Segment seg;
    int interior_sign; // +1: interior is to the left, -1: interior is to the right
};

struct EdgeFragment {
    Point a;
    Point b;
    int interior_sign;
};

std::vector<EdgeInfo> polygon_edges(const Polygon& poly) {
    std::vector<EdgeInfo> edges;
    for (const auto& ring : poly.rings) {
        const int n = static_cast<int>(ring.vertices.size());
        if (n < 2) {
            continue;
        }

        const int interior_sign = (signed_area(ring) >= 0.0) ? 1 : -1;
        for (int i = 0; i < n; ++i) {
            edges.push_back(EdgeInfo{
                Segment{ring.vertices[i], ring.vertices[(i + 1) % n]},
                interior_sign
            });
        }
    }
    return edges;
}

Point midpoint(const EdgeFragment& f) {
    return Point{
        0.5 * (f.a.x + f.b.x),
        0.5 * (f.a.y + f.b.y)
    };
}

void add_unique_point(std::vector<Point>& pts, const Point& p) {
    for (const auto& q : pts) {
        if (point_equal(p, q)) {
            return;
        }
    }
    pts.push_back(p);
}

std::vector<Point> edge_split_points(
    const Segment& e,
    const std::vector<EdgeInfo>& other_edges)
{
    std::vector<Point> pts;
    pts.push_back(e.a);
    pts.push_back(e.b);

    for (const auto& other : other_edges) {
        Point ip{};
        if (segment_intersection_point(e.a, e.b, other.seg.a, other.seg.b, ip)) {
            add_unique_point(pts, ip);
        }

        // Handle collinear-overlap cases by inserting overlapping endpoints.
        if (orientation(e.a, e.b, other.seg.a) == 0 && orientation(e.a, e.b, other.seg.b) == 0) {
            if (on_segment(e.a, e.b, other.seg.a)) {
                add_unique_point(pts, other.seg.a);
            }
            if (on_segment(e.a, e.b, other.seg.b)) {
                add_unique_point(pts, other.seg.b);
            }
        }
    }

    std::sort(pts.begin(), pts.end(),
        [&](const Point& p1, const Point& p2) {
            return segment_parameter(e.a, e.b, p1) < segment_parameter(e.a, e.b, p2);
        });

    std::vector<Point> unique_sorted;
    for (const auto& p : pts) {
        if (unique_sorted.empty() || !point_equal(unique_sorted.back(), p)) {
            unique_sorted.push_back(p);
        }
    }

    return unique_sorted;
}

std::vector<EdgeFragment> split_segment_by_points(
    const EdgeInfo& e,
    const std::vector<Point>& pts)
{
    std::vector<EdgeFragment> out;
    for (size_t i = 0; i + 1 < pts.size(); ++i) {
        if (!point_equal(pts[i], pts[i + 1])) {
            out.push_back(EdgeFragment{pts[i], pts[i + 1], e.interior_sign});
        }
    }
    return out;
}

Point interior_sample_point(const EdgeFragment& frag) {
    const Point m = midpoint(frag);
    const double dx = frag.b.x - frag.a.x;
    const double dy = frag.b.y - frag.a.y;
    const double len = std::sqrt(dx * dx + dy * dy);
    if (len < EPS) {
        return m;
    }

    const double nx_left = -dy / len;
    const double ny_left =  dx / len;
    const double scale = 1e-7 * len + 1e-9;
    const double dir = static_cast<double>(frag.interior_sign);

    return Point{
        m.x + dir * scale * nx_left,
        m.y + dir * scale * ny_left
    };
}

bool fragment_is_kept_for_intersection(
    const EdgeFragment& frag,
    const Polygon& other,
    bool include_shared_boundary)
{
    const Point m = midpoint(frag);
    const bool midpoint_on_other_boundary = point_on_polygon_boundary(m, other);
    const Point interior_probe = interior_sample_point(frag);
    const bool overlap_filled_region = point_in_polygon_with_holes(interior_probe, other);

    if (!overlap_filled_region) {
        return false;
    }

    if (midpoint_on_other_boundary) {
        return include_shared_boundary;
    }

    return true;
}

std::vector<EdgeFragment> kept_intersection_fragments(
    const Polygon& source,
    const Polygon& other,
    bool include_shared_boundary)
{
    const std::vector<EdgeInfo> src_edges = polygon_edges(source);
    const std::vector<EdgeInfo> oth_edges = polygon_edges(other);

    std::vector<EdgeFragment> kept;

    for (const auto& e : src_edges) {
        const std::vector<Point> split_pts = edge_split_points(e.seg, oth_edges);
        const std::vector<EdgeFragment> frags = split_segment_by_points(e, split_pts);

        for (const auto& frag : frags) {
            if (fragment_is_kept_for_intersection(frag, other, include_shared_boundary)) {
                kept.push_back(frag);
            }
        }
    }

    return kept;
}

std::vector<std::vector<Point>> trace_closed_loops(const std::vector<EdgeFragment>& frags)
{
    std::vector<bool> used(frags.size(), false);
    std::vector<std::vector<Point>> loops;

    for (size_t i = 0; i < frags.size(); ++i) {
        if (used[i]) continue;

        std::vector<Point> loop;
        loop.push_back(frags[i].a);
        loop.push_back(frags[i].b);
        used[i] = true;

        Point current = frags[i].b;
        const Point start = frags[i].a;

        bool progressed = true;
        while (progressed && !point_equal(current, start)) {
            progressed = false;

            for (size_t j = 0; j < frags.size(); ++j) {
                if (used[j]) continue;

                if (point_equal(frags[j].a, current)) {
                    loop.push_back(frags[j].b);
                    current = frags[j].b;
                    used[j] = true;
                    progressed = true;
                    break;
                }
            }
        }

        if (loop.size() >= 4 && point_equal(loop.front(), loop.back())) {
            loop.pop_back();
            loops.push_back(loop);
        }
    }

    return loops;
}

} // namespace

Polygon read_csv(const std::string& path) {
    std::ifstream fin(path);
    if (!fin) {
        throw std::runtime_error("Failed to open input file: " + path);
    }

    std::string line;
    if (!std::getline(fin, line)) {
        throw std::runtime_error("Input file is empty");
    }

    std::map<int, std::vector<Point>> grouped;

    while (std::getline(fin, line)) {
        if (line.empty()) continue;

        std::stringstream ss(line);
        std::string token;

        int ring_id = 0;
        int vertex_id = 0;
        double x = 0.0;
        double y = 0.0;

        if (!std::getline(ss, token, ',')) throw std::runtime_error("Bad CSV row");
        ring_id = std::stoi(token);

        if (!std::getline(ss, token, ',')) throw std::runtime_error("Bad CSV row");
        vertex_id = std::stoi(token);
        (void)vertex_id;

        if (!std::getline(ss, token, ',')) throw std::runtime_error("Bad CSV row");
        x = std::stod(token);

        if (!std::getline(ss, token, ',')) throw std::runtime_error("Bad CSV row");
        y = std::stod(token);

        grouped[ring_id].push_back(Point{x, y});
    }

    Polygon poly;
    poly.rings.reserve(grouped.size());
    for (const auto& [id, verts] : grouped) {
        Ring ring{id, verts};
        remove_consecutive_duplicate_points(ring.vertices);
        poly.rings.push_back(std::move(ring));
    }

    return poly;
}

double signed_area(const Ring& ring) {
    return polygon_signed_area(ring.vertices);
}

double total_signed_area(const Polygon& poly) {
    double total = 0.0;
    for (const auto& ring : poly.rings) {
        total += signed_area(ring);
    }
    return total;
}

int total_vertices(const Polygon& poly) {
    int total = 0;
    for (const auto& ring : poly.rings) {
        total += static_cast<int>(ring.vertices.size());
    }
    return total;
}

bool ring_is_simple(const Ring& ring) {
    const auto& v = ring.vertices;
    const int n = static_cast<int>(v.size());
    if (n < 3) return false;

    for (int i = 0; i < n; ++i) {
        if (point_equal(v[i], v[(i + 1) % n])) {
            return false;
        }
    }

    for (int i = 0; i < n; ++i) {
        const Point a = v[i];
        const Point b = v[(i + 1) % n];
        for (int j = i + 1; j < n; ++j) {
            if (edges_are_adjacent(i, j, n)) continue;

            const Point c = v[j];
            const Point d = v[(j + 1) % n];

            if (segments_intersect(a, b, c, d)) {
                return false;
            }
        }
    }
    return true;
}

bool rings_intersect(const Ring& a, const Ring& b) {
    const int na = static_cast<int>(a.vertices.size());
    const int nb = static_cast<int>(b.vertices.size());

    for (int i = 0; i < na; ++i) {
        const Point a1 = a.vertices[i];
        const Point a2 = a.vertices[(i + 1) % na];
        for (int j = 0; j < nb; ++j) {
            const Point b1 = b.vertices[j];
            const Point b2 = b.vertices[(j + 1) % nb];
            if (segments_intersect(a1, a2, b1, b2)) {
                return true;
            }
        }
    }
    return false;
}

bool point_in_ring(const Point& p, const Ring& ring) {
    const auto& v = ring.vertices;
    const int n = static_cast<int>(v.size());
    if (n < 3) return false;

    for (int i = 0; i < n; ++i) {
        const Point a = v[i];
        const Point b = v[(i + 1) % n];
        if (on_segment(a, b, p)) {
            return true;
        }
    }

    bool inside = false;
    for (int i = 0, j = n - 1; i < n; j = i++) {
        const Point& a = v[i];
        const Point& b = v[j];

        const bool crosses =
            ((a.y > p.y) != (b.y > p.y)) &&
            (p.x < (b.x - a.x) * (p.y - a.y) / ((b.y - a.y) + EPS) + a.x);

        if (crosses) {
            inside = !inside;
        }
    }

    return inside;
}

bool point_strictly_in_ring(const Point& p, const Ring& ring) {
    const auto& v = ring.vertices;
    const int n = static_cast<int>(v.size());
    if (n < 3) return false;

    for (int i = 0; i < n; ++i) {
        const Point a = v[i];
        const Point b = v[(i + 1) % n];
        if (on_segment(a, b, p)) {
            return false;
        }
    }

    bool inside = false;
    for (int i = 0, j = n - 1; i < n; j = i++) {
        const Point& a = v[i];
        const Point& b = v[j];

        const bool crosses =
            ((a.y > p.y) != (b.y > p.y)) &&
            (p.x < (b.x - a.x) * (p.y - a.y) / ((b.y - a.y) + EPS) + a.x);

        if (crosses) {
            inside = !inside;
        }
    }

    return inside;
}

bool polygon_topology_valid(const Polygon& poly) {
    if (poly.rings.empty()) return false;

    for (const auto& ring : poly.rings) {
        if (!ring_is_simple(ring)) return false;
    }

    for (std::size_t i = 0; i < poly.rings.size(); ++i) {
        for (std::size_t j = i + 1; j < poly.rings.size(); ++j) {
            if (rings_intersect(poly.rings[i], poly.rings[j])) {
                return false;
            }
        }
    }

    const Ring& exterior = poly.rings[0];
    for (std::size_t i = 1; i < poly.rings.size(); ++i) {
        if (poly.rings[i].vertices.empty()) return false;
        if (!point_strictly_in_ring(poly.rings[i].vertices[0], exterior)) {
            return false;
        }
    }

    for (std::size_t i = 1; i < poly.rings.size(); ++i) {
        for (std::size_t j = i + 1; j < poly.rings.size(); ++j) {
            if (poly.rings[j].vertices.empty()) return false;

            if (point_in_ring(poly.rings[j].vertices[0], poly.rings[i])) {
                return false;
            }
            if (point_in_ring(poly.rings[i].vertices[0], poly.rings[j])) {
                return false;
            }
        }
    }

    return true;
}

bool point_in_polygon_with_holes(const Point& p, const Polygon& poly) {
    if (poly.rings.empty()) return false;

    if (!point_in_ring(p, poly.rings[0])) return false;

    for (std::size_t i = 1; i < poly.rings.size(); ++i) {
        if (point_in_ring(p, poly.rings[i])) {
            return false;
        }
    }
    return true;
}

static Polygon make_single_ring_polygon(const Ring& ring)
{
    Polygon p;
    p.rings.push_back(ring);
    return p;
}

static double simple_ring_intersection_area(const Ring& a, const Ring& b)
{
    const Polygon pa = make_single_ring_polygon(a);
    const Polygon pb = make_single_ring_polygon(b);

    std::vector<EdgeFragment> kept = kept_intersection_fragments(pa, pb, true);
    std::vector<EdgeFragment> kept_other = kept_intersection_fragments(pb, pa, false);
    kept.insert(kept.end(), kept_other.begin(), kept_other.end());

    const std::vector<std::vector<Point>> loops = trace_closed_loops(kept);

    double area = 0.0;
    for (const auto& loop : loops) {
        area += std::fabs(polygon_signed_area(loop));
    }
    return area;
}

double intersection_area_between(const Polygon& a, const Polygon& b)
{
    const std::size_t ring_count = std::min(a.rings.size(), b.rings.size());
    double total = 0.0;
    for (std::size_t i = 0; i < ring_count; ++i) {
        total += simple_ring_intersection_area(a.rings[i], b.rings[i]);
    }
    return total;
}

double total_areal_displacement_between(const Polygon& input, const Polygon& output)
{
    const std::size_t ring_count = std::min(input.rings.size(), output.rings.size());

    double displacement = 0.0;
    for (std::size_t i = 0; i < ring_count; ++i) {
        const double area_in = std::fabs(signed_area(input.rings[i]));
        const double area_out = std::fabs(signed_area(output.rings[i]));
        const double inter = simple_ring_intersection_area(input.rings[i], output.rings[i]);
        displacement += area_in + area_out - 2.0 * inter;
    }

    if (displacement < 0.0 && std::fabs(displacement) < 1e-9) {
        displacement = 0.0;
    }
    return displacement;
}

void normalize_vertex_ids_and_print(const Polygon& output, const Polygon& input, double total_areal_displacement) {
    std::cout << "ring_id,vertex_id,x,y\n";
    std::cout << std::defaultfloat << std::setprecision(10);

    for (std::size_t ring_index = 0; ring_index < output.rings.size(); ++ring_index) {
        const Ring& ring = output.rings[ring_index];
        const int n = static_cast<int>(ring.vertices.size());
        int start = 0;
        if (ring_index < input.rings.size() && !input.rings[ring_index].vertices.empty() && n > 0) {
            const Point& anchor = input.rings[ring_index].vertices[0];
            double best = distance2(ring.vertices[0], anchor);
            for (int i = 1; i < n; ++i) {
                const double d = distance2(ring.vertices[i], anchor);
                if (d < best) {
                    best = d;
                    start = i;
                }
            }
        }

        for (int k = 0; k < n; ++k) {
            const int i = (start + k) % n;
            const Point& p = ring.vertices[i];
            std::cout << ring_index << "," << k << "," << p.x << "," << p.y << "\n";
        }
    }

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Total signed area in input: " << total_signed_area(input) << "\n";
    std::cout << "Total signed area in output: " << total_signed_area(output) << "\n";
    std::cout << "Total areal displacement: " << total_areal_displacement << "\n";
}
