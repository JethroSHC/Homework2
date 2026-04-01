#pragma once

#include <string>
#include <vector>
#include "geometry.h"

/**
 * @struct Ring
 * @brief Represents a single polygon ring (either exterior or interior).
 *
 * A ring is a closed sequence of vertices forming a polygon boundary.
 * - ring_id = 0 -> exterior ring (counterclockwise)
 * - ring_id > 0 -> interior rings / holes (clockwise)
*/
struct Ring {
    int ring_id{};
    std::vector<Point> vertices;
};


/**
 * @struct Polygon
 * @brief Represents a polygon with one exterior ring and zero or more holes.
 *
 * The polygon is composed of multiple rings:
 * - One exterior boundary (ring_id = 0)
 * - Zero or more interior rings (holes)
*/
struct Polygon {
    std::vector<Ring> rings;
};

/**
 * @brief Reads a polygon from a CSV file.
 *
 * The CSV format must follow:
 * ring_id,vertex_id,x,y
 *
 * Vertices are grouped by ring_id.
 *
 * @param path Path to the input CSV file
 * @return Polygon object containing all rings
 * @throws std::runtime_error if file cannot be opened or parsed
*/
Polygon read_csv(const std::string& path);

/**
 * @brief Computes the signed area of a ring.
 *
 * Positive area -> counterclockwise (exterior ring)  
 * Negative area -> clockwise (hole)
 *
 * @param ring Input ring
 * @return Signed area of the ring
*/
double signed_area(const Ring& ring);

/**
 * @brief Computes the total signed area of the polygon.
 *
 * This is the sum of:
 * - exterior ring area (positive)
 * - interior ring areas (negative)
 *
 * @param poly Input polygon
 * @return Total signed area
*/
double total_signed_area(const Polygon& poly);

/**
 * @brief Counts the total number of verts across all rings.
 *
 * @param poly Input polygon
 * @return Total number of verts
*/
int total_vertices(const Polygon& poly);

/**
 * @brief Checks whether a ring is simple (non-self-intersecting).
 *
 * Ensures that no edges intersect except adjacent edges.
 *
 * @param ring Input ring
 * @return true if the ring is simple, false otherwise
*/
bool ring_is_simple(const Ring& ring);

/**
 * @brief Checks whether two rings intersect.
 *
 * Detects if any edge of one ring intersects any edge of another.
 *
 * @param a First ring
 * @param b Second ring
 * @return true if rings intersect, false otherwise
*/
bool rings_intersect(const Ring& a, const Ring& b);

/**
 * @brief Returns true if point p lies inside or on the boundary of ring.
 */
bool point_in_ring(const Point& p, const Ring& ring);

bool point_strictly_in_ring(const Point& p, const Ring& ring);

/**
 * @brief Validates polygon topology:
 * - each ring is simple
 * - no ring intersections
 * - all holes are inside exterior ring
 * - holes do not contain each other
 */
bool polygon_topology_valid(const Polygon& poly);

bool point_in_polygon_with_holes(const Point& p, const Polygon& poly);

double intersection_area_between(const Polygon& a, const Polygon& b);

double total_areal_displacement_between(const Polygon& input, const Polygon& output);

/**
 * @brief Prints the polygon in CSV format with normalized vertex IDs.
 *
 * Output format:
 * ring_id,vertex_id,x,y
 *
 * Also prints:
 * - Total signed area (input and output)
 * - Total areal displacement
 *
 * @param output Simplified polygon
 * @param input Original polygon
 * @param total_areal_displacement Accumulated displacement value
*/
void normalize_vertex_ids_and_print(const Polygon& output, const Polygon& input, double total_areal_displacement);