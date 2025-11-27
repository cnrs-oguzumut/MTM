#pragma once
#ifndef NEIGHBOR_ANALYZER_H
#define NEIGHBOR_ANALYZER_H

#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include "Point2D.h"
#include "src/ap.h"  // ALGLIB header
#include "src/dataanalysis.h"  // ALGLIB KD-tree

struct NeighborChangeInfo {
    bool has_changed;
    int total_nodes_changed;
    int total_connections_added;
    int total_connections_removed;
    std::map<int, std::set<int>> added_neighbors;    // node -> new neighbors
    std::map<int, std::set<int>> removed_neighbors;  // node -> lost neighbors
    std::map<int, int> neighbor_count_change;        // node -> delta in neighbor count
    
    NeighborChangeInfo() 
        : has_changed(false), 
          total_nodes_changed(0),
          total_connections_added(0),
          total_connections_removed(0) {}
};

class NeighborAnalyzer {
public:
    enum class SearchType {
        RADIUS,    // Find all neighbors within a radius
        K_NEAREST  // Find k nearest neighbors
    };
    
    // Constructor
    NeighborAnalyzer(SearchType type = SearchType::K_NEAREST);
    
    // Set search parameters
    void setRadiusSearch(double cutoff_distance);
    void setKNearestSearch(int k_neighbors);
    void setDebugMode(bool enable);
    
    // Build neighbor list from points
    std::map<int, std::set<int>> buildNeighbors(
        const std::vector<Point2D>& points
    ) const;
    
    // Compare two neighbor sets
    static NeighborChangeInfo compareNeighbors(
        const std::map<int, std::set<int>>& neighbors_before,
        const std::map<int, std::set<int>>& neighbors_after
    );
    
    // Print change information
    static void printNeighborChanges(
        const NeighborChangeInfo& info,
        int max_nodes_to_print = 5
    );
    
    // Quick check if neighbors changed (returns only boolean)
    static bool hasNeighborChanged(
        const std::map<int, std::set<int>>& neighbors_before,
        const std::map<int, std::set<int>>& neighbors_after
    );
    
    // Debug: detailed neighbor comparison
    static void debugNeighborComparison(
        const std::vector<Point2D>& points_before,
        const std::vector<Point2D>& points_after,
        const std::map<int, std::set<int>>& neighbors_before,
        const std::map<int, std::set<int>>& neighbors_after,
        int node_to_inspect = -1
    );

    // Dans NeighborAnalyzer.h, ajoutez:
static NeighborChangeInfo compareNeighborsWithTolerance(
    const std::vector<Point2D>& points_before,
    const std::vector<Point2D>& points_after,
    const std::map<int, std::set<int>>& neighbors_before,
    const std::map<int, std::set<int>>& neighbors_after,
    double distance_tolerance = 0.1  // 10% de tolérance
);
    
private:
    SearchType search_type_;
    double cutoff_distance_;
    int k_neighbors_;
    bool debug_mode_;
    
    // Internal methods for different search types
    std::map<int, std::set<int>> buildRadiusNeighbors(
        const std::vector<Point2D>& points
    ) const;
    
    std::map<int, std::set<int>> buildKNearestNeighbors(
        const std::vector<Point2D>& points
    ) const;
};

// Standalone diagnostic function (not part of the class)
void diagnosePointDistribution(const std::vector<Point2D>& points);

#endif // NEIGHBOR_ANALYZER_H