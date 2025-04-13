#include "../include/loading/NanoIndenter.h"
#include <limits>
#include <cmath>

/**
 * Constructor
 */
NanoIndenter::NanoIndenter(double r, const Eigen::Vector2d& initial_pos)
    : radius(r), position(initial_pos), indentation_depth(0.0), initial_position(initial_pos) {
}

/**
 * Initialize indenter position based on the upper layer of atoms
 */
// Add this to your initialization method
void NanoIndenter::initializePosition(const std::vector<Point2D>& atoms) {
    // Step 1: Find domain dimensions
    double min_x = std::numeric_limits<double>::max();
    double max_x = -std::numeric_limits<double>::max();
    double min_y = std::numeric_limits<double>::max();
    double max_y = -std::numeric_limits<double>::max();
    
    for (const auto& atom : atoms) {
        min_x = std::min(min_x, atom.coord.x());
        max_x = std::max(max_x, atom.coord.x());
        min_y = std::min(min_y, atom.coord.y());
        max_y = std::max(max_y, atom.coord.y());
    }
    
    // Step 2: Define a band thickness for the top two layers
    double domain_height = max_y - min_y;
    double top_band_thickness = domain_height * 0.15; // Top 15% of the domain height
    double top_band_lower_bound = max_y - top_band_thickness;
    
    // Step 3: Collect atoms from the top band (which will contain the top two layers)
    std::vector<int> top_band_indices;
    for (int i = 0; i < static_cast<int>(atoms.size()); i++) {
        if (atoms[i].coord.y() >= top_band_lower_bound) {
            top_band_indices.push_back(i);
        }
    }
    
    // Step 4: If we found a reasonable number of atoms, try to identify two distinct layers
    // by clustering their y-coordinates
    top_layer_atom_indices.clear();
    std::vector<int> second_layer_indices;
    
    if (top_band_indices.size() >= 3) {  // Need at least a few atoms to work with
        // Sort the atoms in the band by y-coordinate (highest first)
        std::sort(top_band_indices.begin(), top_band_indices.end(), 
                 [&atoms](int a, int b) { 
                     return atoms[a].coord.y() > atoms[b].coord.y(); 
                 });
        
        // Analyze y-coordinates to find layers
        std::vector<double> y_values;
        for (int idx : top_band_indices) {
            y_values.push_back(atoms[idx].coord.y());
        }
        
        // Look for a significant gap in y-values to identify the layer separation
        double prev_y = y_values[0];
        double largest_gap = 0.0;
        size_t gap_index = 0;
        
        for (size_t i = 1; i < y_values.size(); i++) {
            double gap = prev_y - y_values[i];
            if (gap > largest_gap) {
                largest_gap = gap;
                gap_index = i;
            }
            prev_y = y_values[i];
        }
        
        // If we found a significant gap (compared to the lattice spacing)
        // use it to separate the two layers
        double avg_spacing = top_band_thickness / (y_values.size() - 1);
        if (largest_gap > 1.5 * avg_spacing && gap_index > 0) {
            // First layer: atoms before the gap
            for (size_t i = 0; i < gap_index; i++) {
                top_layer_atom_indices.push_back(top_band_indices[i]);
            }
            
            // Second layer: atoms after the gap
            for (size_t i = gap_index; i < top_band_indices.size(); i++) {
                second_layer_indices.push_back(top_band_indices[i]);
            }
        } else {
            // Couldn't find a clear gap, so divide the band into two equal parts
            size_t mid_point = top_band_indices.size() / 2;
            
            // First layer: top half of atoms
            for (size_t i = 0; i < mid_point; i++) {
                top_layer_atom_indices.push_back(top_band_indices[i]);
            }
            
            // Second layer: bottom half of atoms
            for (size_t i = mid_point; i < top_band_indices.size(); i++) {
                second_layer_indices.push_back(top_band_indices[i]);
            }
        }
    } else {
        // Not enough atoms found, use all atoms in the top band
        for (int idx : top_band_indices) {
            top_layer_atom_indices.push_back(idx);
        }
    }
    
    // Step 5: Combine both layers and find the center
    std::vector<int> combined_indices;
    combined_indices.insert(combined_indices.end(), 
                          top_layer_atom_indices.begin(), 
                          top_layer_atom_indices.end());
    combined_indices.insert(combined_indices.end(), 
                          second_layer_indices.begin(), 
                          second_layer_indices.end());
    
    // Calculate the center position from all atoms in the two layers
    double sum_x = 0.0;
    for (int idx : combined_indices) {
        sum_x += atoms[idx].coord.x();
    }
    
    // Find the x-center of our combined layers
    double center_x = combined_indices.empty() ? (min_x + max_x) / 2.0 : 
                     sum_x / combined_indices.size();
    
    // Step 6: Position the indenter above the center of the combined layers
    position = Eigen::Vector2d(center_x, max_y + radius);
    initial_position = position;
    indentation_depth = 0.0;
    
    std::cout << "Initialized indenter based on top " << combined_indices.size() 
              << " atoms (layer 1: " << top_layer_atom_indices.size() 
              << ", layer 2: " << second_layer_indices.size() << ")" << std::endl;
    std::cout << "Indenter position: (" << position.x() << ", " << position.y() << ")" << std::endl;
}
/**
 * Check if an atom is in contact with the indenter
 */
bool NanoIndenter::isInContact(const Eigen::Vector2d& atom_pos) const {
    // Implement exactly like the working code in main
    double distance = (atom_pos - position).norm();
    
    // Print debug info for troubleshooting
    // std::cout << "Checking atom at (" << atom_pos.x() << ", " << atom_pos.y() 
    //           << "), distance: " << distance << ", radius: " << radius 
    //           << ", in contact: " << (distance <= radius + 1e-10) << std::endl;
              
    return distance <= radius + 1e-10;
}
 /**
 * Project an atom onto the indenter surface
 */
Eigen::Vector2d NanoIndenter::projectOntoSurface(const Eigen::Vector2d& atom_pos) const {
    Eigen::Vector2d direction = (atom_pos - position).normalized();
    return position + direction * radius;
}

/**
 * Get indenter position
 */
const Eigen::Vector2d& NanoIndenter::getPosition() const {
    return position;
}

/**
 * Get indenter radius
 */
double NanoIndenter::getRadius() const {
    return radius;
}

/**
 * Get indices of atoms in contact with the indenter
 */
const std::vector<int>& NanoIndenter::getFixedAtoms() const {
    return fixed_atoms;
}

/**
 * Get current indentation depth
 */
double NanoIndenter::getIndentationDepth() const {
    return indentation_depth;
}

/**
 * Get y-coordinate of the bottom-most point of the indenter
 */
double NanoIndenter::getBottomY() const {
    return position.y() - radius;
}

std::vector<int> NanoIndenter::findTopLayerAtoms(const std::vector<Point2D>& atoms) const {
    // Find the maximum y-coordinate
    double max_y = -std::numeric_limits<double>::infinity();
    for (const auto& atom : atoms) {
        if (atom.coord.y() > max_y) {
            max_y = atom.coord.y();
        }
    }
    
    // Collect indices of atoms in the top layer
    std::vector<int> top_layer_indices;
    const double tolerance = 1e-6;
    
    for (size_t i = 0; i < atoms.size(); i++) {
        if (std::abs(atoms[i].coord.y() - max_y) < tolerance) {
            top_layer_indices.push_back(i);
        }
    }
    
    return top_layer_indices;
}

// Then modify setIndentationDepth to use the stored indices
void NanoIndenter::setIndentationDepth(double depth, const std::vector<Point2D>& atoms) {
    indentation_depth = depth;
    
    // Update indenter position
    position = initial_position;
    position.y() -= indentation_depth;
    
    // Clear previous contacts
    fixed_atoms.clear();
    
    // Check for contacts only with top layer atoms (using stored indices)
    for (int idx : top_layer_atom_indices) {
        if (isInContact(atoms[idx].coord)) {
            fixed_atoms.push_back(idx);
        }
    }
}
void NanoIndenter::updateFixedAtomPositions(std::vector<Point2D>& atoms) {
    // For each atom in contact with the indenter
    for (int idx : fixed_atoms) {
        // Project the atom onto the indenter surface
        atoms[idx].coord = projectOntoSurface(atoms[idx].coord);
    }
}

std::vector<int> NanoIndenter::getNeighboringAtomIndices(const std::vector<Point2D>& atoms, 
    const std::vector<int>& contact_atoms,
    double horizontal_tolerance ,
    double vertical_tolerance ) {
std::vector<int> neighboring_atoms;  // To store the indices of neighboring atoms

// Loop over all the contact atoms to find neighboring atoms below them
for (int contact_idx : contact_atoms) {
const auto& contact_atom = atoms[contact_idx];

// Loop over all atoms and find those that are nearby and below the contact atom
for (size_t i = 0; i < atoms.size(); ++i) {
if (i == contact_idx) continue;  // Skip the contact atom itself

const auto& atom = atoms[i];
double dx = atom.coord.x() - contact_atom.coord.x();
double dy = atom.coord.y() - contact_atom.coord.y();

// Check if the atom is below the contact atom and within tolerance
if (std::abs(dx) < horizontal_tolerance && dy < 0 && std::abs(dy) < vertical_tolerance) {
neighboring_atoms.push_back(i);
}
}
}

return neighboring_atoms;
}

Eigen::Vector2d NanoIndenter::projectOntoSurface_NN(const Eigen::Vector2d& atom_pos) const {
    // Calculate the direction from the atom position to the indenter surface
    Eigen::Vector2d direction = (atom_pos - position).normalized();
    
    // Apply the projection with the original radius, with an additional offset to avoid overlap
    return position + direction * (radius + 0.2); // dR is a small offset to avoid superposition
}
