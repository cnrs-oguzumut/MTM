#include "../include/mesh/MeshGenerator.h"
#include "../include/output/triangulation_to_vtk.h"
#include <sys/stat.h>

void write_triangulation_to_vtk(
    const std::vector<Triangle>& triangles,
    const std::vector<Point2D>& points,
    int iteration_number,
    const std::string& prefix
) {
    // Create vtk directory if it doesn't exist
    const std::string vtk_dir = "vtk";
    struct stat info;
    if (stat(vtk_dir.c_str(), &info) != 0) {
        // Directory doesn't exist, create it
        int result = mkdir(vtk_dir.c_str(), 0777);
        if (result != 0) {
            std::cerr << "Error: Could not create directory " << vtk_dir << std::endl;
            return;
        }
    }
    
    // Create filename with directory
    std::string filename = vtk_dir + "/" + prefix + "_" + std::to_string(iteration_number) + ".vtk";
    std::ofstream vtk_file(filename);
    
    if (!vtk_file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing!" << std::endl;
        return;
    }
    
    // Rest of the function remains the same
    // Extract unique node indices from the triangles
    std::cout << "Creating VTK file: extracting unique nodes... " << std::flush;
    std::unordered_set<int> unique_indices;
    for (const auto& tri : triangles) {
        unique_indices.insert(tri.vertex_indices[0]);
        unique_indices.insert(tri.vertex_indices[1]);
        unique_indices.insert(tri.vertex_indices[2]);
    }
    std::cout << unique_indices.size() << " unique nodes found" << std::endl;
    
    // Create a mapping from old indices to new indices
    std::vector<int> indices_vector(unique_indices.begin(), unique_indices.end());
    std::sort(indices_vector.begin(), indices_vector.end());
    std::unordered_map<int, int> index_map;
    for (size_t i = 0; i < indices_vector.size(); ++i) {
        index_map[indices_vector[i]] = i;
    }
    
    // Begin writing the VTK file
    vtk_file << "# vtk DataFile Version 1.0\n";
    vtk_file << "2D Unstructured Grid of Linear Triangles\n";
    vtk_file << "ASCII\n";
    vtk_file << " \n";
    vtk_file << "DATASET UNSTRUCTURED_GRID\n";
    vtk_file << "POINTS " << unique_indices.size() << " float\n";
    
    // Write only the unique points
    std::cout << "Writing points to VTK... " << std::flush;
    for (int idx : indices_vector) {
        vtk_file << points[idx].coord(0) << " " << points[idx].coord(1) << " 0\n";
    }
    std::cout << "done" << std::endl;
    
    // Write cells (triangles)
    vtk_file << "CELLS " << triangles.size() << " " << triangles.size() * 4 << "\n";
    for (const auto& tri : triangles) {
        vtk_file << "3 "
                << index_map[tri.vertex_indices[0]] << " "
                << index_map[tri.vertex_indices[1]] << " "
                << index_map[tri.vertex_indices[2]] << "\n";
    }
    
    // Write cell types
    std::cout << "Finalizing VTK file... " << std::flush;
    vtk_file << "CELL_TYPES " << triangles.size() << "\n";
    for (int id = 0; id < triangles.size(); id++) {
        vtk_file << "5\n"; // Type 5 is a triangle
    }
    
    vtk_file.close();
    std::cout << "done - VTK file created successfully: " << filename << std::endl;
}