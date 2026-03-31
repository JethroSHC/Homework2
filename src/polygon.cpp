#include "polygon.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <map>

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
    for (const auto& [id, verts] : grouped) {
        poly.rings.push_back(Ring{id, verts});
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

static bool edges_are_adjacent(int i, int j, int n) {
    if (i == j) return true;
    if ((i + 1) % n == j) return true;
    if ((j + 1) % n == i) return true;
    return false;
}

bool ring_is_simple(const Ring& ring) {
    const auto& v = ring.vertices;
    const int n = static_cast<int>(v.size());
    if (n < 3) return false;

    for (int i = 0; i < n; ++i) {
        Point a = v[i];
        Point b = v[(i + 1) % n];
        for (int j = i + 1; j < n; ++j) {
            if (edges_are_adjacent(i, j, n)) continue;

            Point c = v[j];
            Point d = v[(j + 1) % n];

            // allow shared endpoint only for adjacency, already skipped above
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
        Point a1 = a.vertices[i];
        Point a2 = a.vertices[(i + 1) % na];
        for (int j = 0; j < nb; ++j) {
            Point b1 = b.vertices[j];
            Point b2 = b.vertices[(j + 1) % nb];
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

    // Boundary counts as inside
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

bool polygon_topology_valid(const Polygon& poly) {
    if (poly.rings.empty()) return false;

    // 1. Every ring must be simple
    for (const auto& ring : poly.rings) {
        if (!ring_is_simple(ring)) return false;
    }

    // 2. No pair of rings may intersect
    for (size_t i = 0; i < poly.rings.size(); ++i) {
        for (size_t j = i + 1; j < poly.rings.size(); ++j) {
            if (rings_intersect(poly.rings[i], poly.rings[j])) {
                return false;
            }
        }
    }

    // 3. Every interior ring must lie inside the exterior ring
    const Ring& exterior = poly.rings[0];
    for (size_t i = 1; i < poly.rings.size(); ++i) {
        if (poly.rings[i].vertices.empty()) return false;
        if (!point_in_ring(poly.rings[i].vertices[0], exterior)) {
            return false;
        }
    }

    // 4. Holes must not contain each other
    for (size_t i = 1; i < poly.rings.size(); ++i) {
        for (size_t j = i + 1; j < poly.rings.size(); ++j) {
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


void normalize_vertex_ids_and_print(const Polygon& output, const Polygon& input, double total_areal_displacement) {
    std::cout << "ring_id,vertex_id,x,y\n";
    for (const auto& ring : output.rings) {
        for (int i = 0; i < static_cast<int>(ring.vertices.size()); ++i) {
            const Point& p = ring.vertices[i];
            std::cout << ring.ring_id << "," << i << "," << p.x << "," << p.y << "\n";
        }
    }

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Total signed area in input: " << total_signed_area(input) << "\n";
    std::cout << "Total signed area in output: " << total_signed_area(output) << "\n";
    std::cout << "Total areal displacement: " << total_areal_displacement << "\n";
}