#ifndef NANOINDENTER_H
#define NANOINDENTER_H

#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <algorithm>

#include "../include/geometry/Point2D.h"


/**
 * @brief Class for simulating nanoindentation experiments
 * 
 * This class models a spherical indenter that can be positioned above a material
 * surface and pressed into it. It handles contact detection, atom fixation,
 * and surface projection for atoms in contact with the indenter.
 */
class NanoIndenter {
private:
    double radius;                           ///< Radius of the spherical indenter
    Eigen::Vector2d position;                ///< Current position of indenter center
    Eigen::Vector2d initial_position;        ///< Initial position of indenter center
    double indentation_depth;                ///< Current indentation depth
    std::vector<int> fixed_atoms;            ///< Indices of atoms in contact with indenter
    std::vector<int> top_layer_atom_indices; ///< Indices of atoms in the top layer

public:
    /**
     * @brief Constructor
     * @param r Radius of the spherical indenter
     * @param initial_pos Initial position of the indenter center
     */
    NanoIndenter(double r, const Eigen::Vector2d& initial_pos);

    /**
     * @brief Initialize indenter position based on the upper layer of atoms
     * 
     * Automatically positions the indenter above the center of the top two
     * atomic layers. Uses clustering analysis to identify distinct layers.
     * 
     * @param atoms Vector of all atoms in the system
     */
    void initializePosition(const std::vector<Point2D>& atoms);

    /**
     * @brief Check if an atom is in contact with the indenter
     * @param atom_pos Position of the atom to check
     * @return true if the atom is within the indenter radius
     */
    bool isInContact(const Eigen::Vector2d& atom_pos) const;

    /**
     * @brief Project an atom onto the indenter surface
     * @param atom_pos Current position of the atom
     * @return New position on the indenter surface
     */
    Eigen::Vector2d projectOntoSurface(const Eigen::Vector2d& atom_pos) const;

    /**
     * @brief Project atom onto surface with additional offset to prevent overlap
     * @param atom_pos Current position of the atom
     * @return New position with small offset from indenter surface
     */
    Eigen::Vector2d projectOntoSurface_NN(const Eigen::Vector2d& atom_pos) const;

    /**
     * @brief Set the indentation depth and update contacts
     * @param depth New indentation depth
     * @param atoms Vector of all atoms to check for contacts
     */
    void setIndentationDepth(double depth, const std::vector<Point2D>& atoms);

    /**
     * @brief Update positions of atoms in contact with the indenter
     * @param atoms Vector of atoms to update (modified in place)
     */
    void updateFixedAtomPositions(std::vector<Point2D>& atoms);

    /**
     * @brief Find atoms in the topmost layer
     * @param atoms Vector of all atoms
     * @return Vector of indices of atoms in the top layer
     */
    std::vector<int> findTopLayerAtoms(const std::vector<Point2D>& atoms) const;

    /**
     * @brief Find neighboring atoms below contact atoms
     * @param atoms Vector of all atoms
     * @param contact_atoms Indices of atoms in contact with indenter
     * @param horizontal_tolerance Maximum horizontal distance to consider as neighbor
     * @param vertical_tolerance Maximum vertical distance to consider as neighbor
     * @return Vector of indices of neighboring atoms
     */
    std::vector<int> getNeighboringAtomIndices(
        const std::vector<Point2D>& atoms,
        const std::vector<int>& contact_atoms,
        double horizontal_tolerance = 2.0,
        double vertical_tolerance = 2.0);

    // Getters
    /**
     * @brief Get current indenter position
     * @return Reference to the position vector
     */
    const Eigen::Vector2d& getPosition() const;

    /**
     * @brief Get indenter radius
     * @return Radius value
     */
    double getRadius() const;

    /**
     * @brief Get indices of atoms currently in contact with indenter
     * @return Reference to vector of fixed atom indices
     */
    const std::vector<int>& getFixedAtoms() const;

    /**
     * @brief Get current indentation depth
     * @return Current depth value
     */
    double getIndentationDepth() const;

    /**
     * @brief Get y-coordinate of the bottommost point of the indenter
     * @return Y-coordinate of indenter bottom
     */
    double getBottomY() const;
};

#endif // NANOINDENTER_H