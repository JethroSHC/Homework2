#include "geometry.h"
#include <limits>

double cross(const Point& a, const Point& b) {
    return a.x * b.y - a.y * b.x;
}

double cross(const Point& a, const Point& b, const Point& c) {
    return cross(Point{b.x - a.x, b.y - a.y}, Point{c.x - a.x, c.y - a.y});
}

double dot(const Point& a, const Point& b) {
    return a.x * b.x + a.y * b.y;
}

double distance2(const Point& a, const Point& b) {
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    return dx * dx + dy * dy;
}

bool almost_equal(double a, double b, double eps) {
    return std::fabs(a - b) <= eps;
}

bool point_equal(const Point& a, const Point& b, double eps) {
    return almost_equal(a.x, b.x, eps) && almost_equal(a.y, b.y, eps);
}

int orientation(const Point& a, const Point& b, const Point& c) {
    const double val = cross(a, b, c);
    if (val > EPS) return 1;
    if (val < -EPS) return -1;
    return 0;
}

bool on_segment(const Point& a, const Point& b, const Point& p) {
    if (orientation(a, b, p) != 0) return false;
    return (std::min(a.x, b.x) - EPS <= p.x && p.x <= std::max(a.x, b.x) + EPS &&
            std::min(a.y, b.y) - EPS <= p.y && p.y <= std::max(a.y, b.y) + EPS);
}

bool segments_intersect(const Point& a, const Point& b, const Point& c, const Point& d) {
    const int o1 = orientation(a, b, c);
    const int o2 = orientation(a, b, d);
    const int o3 = orientation(c, d, a);
    const int o4 = orientation(c, d, b);

    if (o1 != o2 && o3 != o4) return true;

    if (o1 == 0 && on_segment(a, b, c)) return true;
    if (o2 == 0 && on_segment(a, b, d)) return true;
    if (o3 == 0 && on_segment(c, d, a)) return true;
    if (o4 == 0 && on_segment(c, d, b)) return true;

    return false;
}

double polygon_signed_area(const std::vector<Point>& pts) {
    const int n = static_cast<int>(pts.size());
    if (n < 3) return 0.0;

    double s = 0.0;
    for (int i = 0; i < n; ++i) {
        const Point& p = pts[i];
        const Point& q = pts[(i + 1) % n];
        s += cross(p, q);
    }
    return 0.5 * s;
}

double triangle_signed_area2(const Point& a, const Point& b, const Point& c) {
    return cross(a, b) + cross(b, c) + cross(c, a);
}

static bool line_intersection_general(
    double a1, double b1, double c1,
    double a2, double b2, double c2,
    Point& out)
{
    const double det = a1 * b2 - a2 * b1;
    if (std::fabs(det) < EPS) {
        return false;
    }

    out.x = (b1 * c2 - b2 * c1) / det;
    out.y = (a2 * c1 - a1 * c2) / det;
    return true;
}

static bool line_through_points(const Point& P, const Point& Q, double& a, double& b, double& c) {
    a = Q.y - P.y;
    b = P.x - Q.x;
    c = -(a * P.x + b * P.y);
    return !(std::fabs(a) < EPS && std::fabs(b) < EPS);
}

ApscPlacement apsc_placement_and_displacement(
    const Point& A,
    const Point& B,
    const Point& C,
    const Point& D)
{
    // a*xE + b*yE + c = 0
    const double a = D.y - A.y;
    const double b = A.x - D.x;
    const double c =
        -B.y * A.x +
        (A.y - C.y) * B.x +
        (B.y - D.y) * C.x +
        C.y * D.x;

    auto side_of_directed_line = [](const Point& P, const Point& U, const Point& V) -> int {
        const double s = cross(U, V, P);
        if (s > EPS) return 1;
        if (s < -EPS) return -1;
        return 0;
    };

    auto dist_num_to_line = [](const Point& P, const Point& U, const Point& V) -> double {
        // numerator of perpendicular distance; denominator is constant for both B and C
        return std::fabs(cross(U, V, P));
    };

    // Singular case from paper Fig. 6(b): B,C,D collinear => optimal E = B, displacement 0
    if (std::fabs(cross(B, C, D)) < EPS) {
        return {B, 0.0, true};
    }

    // Build candidate intersections
    Point Eab{};
    bool haveEab = false;
    {
        double a2, b2, c2;
        if (line_through_points(A, B, a2, b2, c2)) {
            haveEab = line_intersection_general(a, b, c, a2, b2, c2, Eab);
        }
    }

    Point Ecd{};
    bool haveEcd = false;
    {
        double a2, b2, c2;
        if (line_through_points(C, D, a2, b2, c2)) {
            haveEcd = line_intersection_general(a, b, c, a2, b2, c2, Ecd);
        }
    }

    // Pick any point on Ē to determine which side of AD it lies on
    Point ElinePoint{};
    bool haveElinePoint = false;
    if (std::fabs(b) > EPS) {
        ElinePoint = Point{0.0, -c / b};
        haveElinePoint = true;
    } else if (std::fabs(a) > EPS) {
        ElinePoint = Point{-c / a, 0.0};
        haveElinePoint = true;
    }

    const int sideB = side_of_directed_line(B, A, D);
    const int sideC = side_of_directed_line(C, A, D);
    const int sideEline = haveElinePoint ? side_of_directed_line(ElinePoint, A, D) : 0;

    Point chosenE{};
    bool chosen = false;

    // Singular case from paper Fig. 6(a): Ē coincident with AD
    if (sideEline == 0) {
        // any point on AD is optimal; choose A
        return {A, 0.0, true};
    }

    // Paper placement pseudocode, page 12
    if (sideB == sideC) {
        if (dist_num_to_line(B, A, D) > dist_num_to_line(C, A, D)) {
            if (haveEab) {
                chosenE = Eab;
                chosen = true;
            }
        } else {
            if (haveEcd) {
                chosenE = Ecd;
                chosen = true;
            }
        }
    } else {
        if (sideB == sideEline) {
            if (haveEab) {
                chosenE = Eab;
                chosen = true;
            }
        } else {
            if (haveEcd) {
                chosenE = Ecd;
                chosen = true;
            }
        }
    }

    // Robust fallback if the preferred line intersection was unavailable
    if (!chosen) {
        if (haveEab) {
            chosenE = Eab;
            chosen = true;
        } else if (haveEcd) {
            chosenE = Ecd;
            chosen = true;
        } else {
            return {A, 0.0, false};
        }
    }

    // Compute displacement for the chosen case
    double disp = 0.0;
    if (chosen && haveEab && point_equal(chosenE, Eab)) {
        // displaced region: E-B-C-D
        disp =
            0.5 * std::fabs(triangle_signed_area2(Eab, B, C)) +
            0.5 * std::fabs(triangle_signed_area2(Eab, C, D));
    } else {
        // displaced region: A-B-C-E
        disp =
            0.5 * std::fabs(triangle_signed_area2(A, B, C)) +
            0.5 * std::fabs(triangle_signed_area2(A, C, chosenE));
    }

    return {chosenE, disp, true};
}

Point area_preserving_point_apsc(
    const Point& A,
    const Point& B,
    const Point& C,
    const Point& D)
{
    return apsc_placement_and_displacement(A, B, C, D).E;
}

double local_displacement_proxy(const Point& A, const Point& B, const Point& C, const Point& D, const Point& E) {
    const std::vector<Point> poly{A, B, C, D, E};
    return std::fabs(polygon_signed_area(poly));
}