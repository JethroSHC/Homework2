#pragma once

#include <vector>
#include <cmath>

/***
 * @struct Point
 * @brief Represents a 2D point.
*/
struct Point {
    double x{};
    double y{};
};

/**
 * @struct Segment
 * @brief Represents a line segment between two points.
*/
struct Segment {
    Point a;
    Point b;
};

/**
 * @brief For floating pt comparison
*/
constexpr double EPS = 1e-9;

/**
 * @brief Computes cross product of two vect
 *
 * @param a First vect
 * @param b Second vect
 * @return Cross product (a.x * b.y - a.y * b.x)
*/
double cross(const Point& a, const Point& b);

/**
 * @brief Computes cross product of vect AB and AC.
 *
 * @param a Point A
 * @param b Point B
 * @param c Point C
 * @return Cross product of AB × AC
*/
double cross(const Point& a, const Point& b, const Point& c);

/**
 * @brief Computes dot product of two vect.
*/
double dot(const Point& a, const Point& b);
double distance2(const Point& a, const Point& b);

bool almost_equal(double a, double b, double eps = EPS);
bool point_equal(const Point& a, const Point& b, double eps = EPS);

/**
 * @brief Determines orientation of three points.
 *
 * @return
 *  1 -> counterclockwise  
 * -1 -> clockwise  
 *  0 -> collinear
*/
int orientation(const Point& a, const Point& b, const Point& c);

/**
 * @brief Checks if pt p lies on segment AB.
*/
bool on_segment(const Point& a, const Point& b, const Point& p);
bool segments_intersect(const Point& a, const Point& b, const Point& c, const Point& d);

bool segment_intersection_point(
    const Point& a, const Point& b,
    const Point& c, const Point& d,
    Point& out);

double segment_parameter(const Point& a, const Point& b, const Point& p);

/**
 * @brief Computes signed area of a polygon (shoelace formula).
 *
 * @param pts List of vertices (closed implicitly)
 * @return Signed area
*/
double polygon_signed_area(const std::vector<Point>& pts);

/**
 * @brief Computes twice the signed area of triangle ABC.
*/
double triangle_signed_area2(const Point& a, const Point& b, const Point& c);

struct ApscPlacement {
    Point E;
    double displacement;
    bool valid;
};

ApscPlacement apsc_placement_and_displacement(
    const Point& A,
    const Point& B,
    const Point& C,
    const Point& D
);

Point area_preserving_point_apsc(
    const Point& A,
    const Point& B,
    const Point& C,
    const Point& D);

// baseline local displacement proxy
double local_displacement_proxy(const Point& A, const Point& B, const Point& C, const Point& D, const Point& E);