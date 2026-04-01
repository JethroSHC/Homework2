#include <iostream>
#include <string>
#include <stdexcept>

#include "polygon.h"
#include "simplify.h"

int main(int argc, char* argv[]) {
    // Validate cmd-line arguments
    if (argc != 3) {
        std::cerr << "Usage: ./simplify <input_file.csv> <target_vertices>\n";
        return 1;
    }

    try {
         // Parse args
        const std::string input_file = argv[1];
        const int target_vertices = std::stoi(argv[2]);

        // Read input polygon
        Polygon input = read_csv(input_file);

        double accumulated_local_displacement = 0.0;
        Polygon output = simplify_polygon(input, target_vertices, accumulated_local_displacement);

        double total_areal_displacement =
            total_areal_displacement_between(input, output);


        //print in format
        normalize_vertex_ids_and_print(output, input, total_areal_displacement);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
