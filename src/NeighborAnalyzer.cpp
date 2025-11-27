#include "../include/geometry/NeighborAnalyzer.h"

NeighborAnalyzer::NeighborAnalyzer(SearchType type)
    : search_type_(type), cutoff_distance_(1.5), k_neighbors_(6), debug_mode_(false) {
}

void NeighborAnalyzer::setRadiusSearch(double cutoff_distance) {
    search_type_ = SearchType::RADIUS;
    cutoff_distance_ = cutoff_distance;
}

void NeighborAnalyzer::setKNearestSearch(int k_neighbors) {
    search_type_ = SearchType::K_NEAREST;
    k_neighbors_ = k_neighbors;
}

void NeighborAnalyzer::setDebugMode(bool enable) {
    debug_mode_ = enable;
}

std::map<int, std::set<int>> NeighborAnalyzer::buildNeighbors(
    const std::vector<Point2D>& points
) const {
    if (search_type_ == SearchType::RADIUS) {
        return buildRadiusNeighbors(points);
    } else {
        return buildKNearestNeighbors(points);
    }
}

std::map<int, std::set<int>> NeighborAnalyzer::buildRadiusNeighbors(
    const std::vector<Point2D>& points
) const {
    std::map<int, std::set<int>> neighbors;
    const int n = points.size();
    
    if (n == 0) return neighbors;
    
    // Build ALGLIB KD-tree with TAGS
    alglib::real_2d_array xy = Point2D::to_alglib_array(points);
    alglib::integer_1d_array tags;
    tags.setlength(n);
    for (int i = 0; i < n; ++i) {
        tags[i] = i;
    }
    
    alglib::kdtree kdt;
    alglib::kdtreebuildtagged(xy, tags, n, 2, 0, 2, kdt);
    
    alglib::kdtreerequestbuffer buf;
    alglib::kdtreecreaterequestbuffer(kdt, buf);
    
    for (int i = 0; i < n; ++i) {
        alglib::real_1d_array query_point = points[i].to_alglib_point();
        
        alglib::integer_1d_array result_tags;
        
        alglib::ae_int_t k_found = alglib::kdtreetsqueryaknn(
            kdt, buf, query_point, cutoff_distance_, 0, true
        );
        alglib::kdtreetsqueryresultstags(kdt, buf, result_tags);
        
        for (int j = 0; j < result_tags.length(); ++j) {
            int neighbor_idx = result_tags[j];
            if (neighbor_idx != i) {
                neighbors[i].insert(neighbor_idx);
            }
        }
    }
    
    return neighbors;
}

std::map<int, std::set<int>> NeighborAnalyzer::buildKNearestNeighbors(
    const std::vector<Point2D>& points
) const {
    std::map<int, std::set<int>> neighbors;
    const int n = points.size();
    
    if (n == 0) return neighbors;
    
    if (debug_mode_) {
        std::cout << "\n=== KD-TREE K-NEAREST SEARCH ===" << std::endl;
        std::cout << "Building KD-tree with " << n << " points" << std::endl;
        std::cout << "Searching for k=" << k_neighbors_ << " nearest neighbors" << std::endl;
    }
    
    // Build ALGLIB KD-tree with TAGS
    alglib::real_2d_array xy = Point2D::to_alglib_array(points);
    alglib::integer_1d_array tags;
    tags.setlength(n);
    for (int i = 0; i < n; ++i) {
        tags[i] = i;
    }
    
    alglib::kdtree kdt;
    alglib::kdtreebuildtagged(xy, tags, n, 2, 0, 2, kdt);
    
    alglib::kdtreerequestbuffer buf;
    alglib::kdtreecreaterequestbuffer(kdt, buf);
    
    for (int i = 0; i < n; ++i) {
        alglib::real_1d_array query_point = points[i].to_alglib_point();
        
        alglib::integer_1d_array result_tags;
        alglib::real_2d_array result_xy;
        
        alglib::ae_int_t k_found = alglib::kdtreetsqueryknn(
            kdt, buf, query_point, k_neighbors_ + 1, true
        );
        alglib::kdtreetsqueryresultstags(kdt, buf, result_tags);
        alglib::kdtreetsqueryresultsxy(kdt, buf, result_xy);
        
        if (debug_mode_ && i < 5) {
            std::cout << "\nNode " << i << " at (" << points[i].x() << ", " << points[i].y() << "):" << std::endl;
            std::cout << "  Requested: " << (k_neighbors_ + 1) << " neighbors" << std::endl;
            std::cout << "  Found: " << k_found << " neighbors" << std::endl;
            std::cout << "  result_tags.length(): " << result_tags.length() << std::endl;
            
            for (int j = 0; j < result_tags.length(); ++j) {
                int tag = result_tags[j];
                double rx = result_xy[j][0];
                double ry = result_xy[j][1];
                double dx = rx - points[i].x();
                double dy = ry - points[i].y();
                double dist = std::sqrt(dx*dx + dy*dy);
                std::cout << "    " << j << ". tag=" << tag << " pos=(" << rx << ", " << ry 
                          << ") dist=" << dist << std::endl;
            }
        }
        
        for (int j = 0; j < result_tags.length(); ++j) {
            int neighbor_idx = result_tags[j];
            if (neighbor_idx != i) {
                neighbors[i].insert(neighbor_idx);
            }
        }
    }
    
    if (debug_mode_) {
        std::cout << "\nTotal nodes with neighbors: " << neighbors.size() << std::endl;
        std::cout << "================================\n" << std::endl;
    }
    
    return neighbors;
}

NeighborChangeInfo NeighborAnalyzer::compareNeighbors(
    const std::map<int, std::set<int>>& neighbors_before,
    const std::map<int, std::set<int>>& neighbors_after
) {
    NeighborChangeInfo info;
    
    std::set<int> all_nodes;
    for (const auto& [node, _] : neighbors_before) {
        all_nodes.insert(node);
    }
    for (const auto& [node, _] : neighbors_after) {
        all_nodes.insert(node);
    }
    
    for (int node : all_nodes) {
        auto it_before = neighbors_before.find(node);
        auto it_after = neighbors_after.find(node);
        
        std::set<int> before_set = (it_before != neighbors_before.end()) 
                                    ? it_before->second 
                                    : std::set<int>();
        std::set<int> after_set = (it_after != neighbors_after.end()) 
                                   ? it_after->second 
                                   : std::set<int>();
        
        std::set<int> added;
        std::set_difference(
            after_set.begin(), after_set.end(),
            before_set.begin(), before_set.end(),
            std::inserter(added, added.begin())
        );
        
        std::set<int> removed;
        std::set_difference(
            before_set.begin(), before_set.end(),
            after_set.begin(), after_set.end(),
            std::inserter(removed, removed.begin())
        );
        
        if (!added.empty() || !removed.empty()) {
            info.has_changed = true;
            info.total_nodes_changed++;
            info.total_connections_added += added.size();
            info.total_connections_removed += removed.size();
            
            if (!added.empty()) {
                info.added_neighbors[node] = added;
            }
            if (!removed.empty()) {
                info.removed_neighbors[node] = removed;
            }
            
            int count_change = static_cast<int>(after_set.size()) - 
                              static_cast<int>(before_set.size());
            info.neighbor_count_change[node] = count_change;
        }
    }
    
    return info;
}

void NeighborAnalyzer::printNeighborChanges(
    const NeighborChangeInfo& info,
    int max_nodes_to_print
) {
    if (!info.has_changed) {
        std::cout << "No neighbor connectivity changes detected." << std::endl;
        return;
    }
    
    std::cout << "Neighbor connectivity CHANGED:" << std::endl;
    std::cout << "  Nodes affected: " << info.total_nodes_changed << std::endl;
    std::cout << "  Total connections added: " << info.total_connections_added << std::endl;
    std::cout << "  Total connections removed: " << info.total_connections_removed << std::endl;
    
    std::vector<std::pair<int, int>> changes;
    for (const auto& [node, delta] : info.neighbor_count_change) {
        changes.emplace_back(node, std::abs(delta));
    }
    std::sort(changes.begin(), changes.end(), 
              [](const auto& a, const auto& b) { return a.second > b.second; });
    
    int n_print = std::min(max_nodes_to_print, static_cast<int>(changes.size()));
    if (n_print > 0) {
        std::cout << "  Most affected nodes (top " << n_print << "):" << std::endl;
        for (int i = 0; i < n_print; ++i) {
            int node = changes[i].first;
            int delta = info.neighbor_count_change.at(node);
            std::cout << "    Node " << node << ": " 
                      << (delta > 0 ? "+" : "") << delta 
                      << " neighbors";
            
            if (info.added_neighbors.count(node)) {
                std::cout << " (added: " << info.added_neighbors.at(node).size() << ")";
            }
            if (info.removed_neighbors.count(node)) {
                std::cout << " (removed: " << info.removed_neighbors.at(node).size() << ")";
            }
            std::cout << std::endl;
        }
    }
}

bool NeighborAnalyzer::hasNeighborChanged(
    const std::map<int, std::set<int>>& neighbors_before,
    const std::map<int, std::set<int>>& neighbors_after
) {
    if (neighbors_before.size() != neighbors_after.size()) {
        return true;
    }
    
    for (const auto& [node, old_neighs] : neighbors_before) {
        auto it = neighbors_after.find(node);
        if (it == neighbors_after.end()) {
            return true;
        }
        if (old_neighs != it->second) {
            return true;
        }
    }
    
    return false;
}

void NeighborAnalyzer::debugNeighborComparison(
    const std::vector<Point2D>& points_before,
    const std::vector<Point2D>& points_after,
    const std::map<int, std::set<int>>& neighbors_before,
    const std::map<int, std::set<int>>& neighbors_after,
    int node_to_inspect
) {
    std::cout << "\n=== NEIGHBOR DEBUG INFO ===" << std::endl;
    std::cout << "Number of points: " << points_before.size() << std::endl;
    std::cout << "Nodes with neighbors (before): " << neighbors_before.size() << std::endl;
    std::cout << "Nodes with neighbors (after): " << neighbors_after.size() << std::endl;
    
    double max_displacement = 0.0;
    double avg_displacement = 0.0;
    int n_moved = 0;
    
    for (size_t i = 0; i < points_before.size() && i < points_after.size(); ++i) {
        double dx = points_after[i].x() - points_before[i].x();
        double dy = points_after[i].y() - points_before[i].y();
        double disp = std::sqrt(dx*dx + dy*dy);
        
        if (disp > 1e-10) {
            n_moved++;
            avg_displacement += disp;
            max_displacement = std::max(max_displacement, disp);
        }
    }
    
    if (n_moved > 0) {
        avg_displacement /= n_moved;
    }
    
    std::cout << "Points moved: " << n_moved << " / " << points_before.size() << std::endl;
    std::cout << "Max displacement: " << max_displacement << std::endl;
    std::cout << "Avg displacement: " << avg_displacement << std::endl;
    
    std::vector<int> nodes_to_check;
    if (node_to_inspect >= 0) {
        nodes_to_check.push_back(node_to_inspect);
    } else {
        for (int i = 0; i < std::min(5, static_cast<int>(points_before.size())); ++i) {
            nodes_to_check.push_back(i);
        }
    }
    
    std::cout << "\n=== DETAILED NODE INSPECTION ===" << std::endl;
    for (int node : nodes_to_check) {
        auto it_before = neighbors_before.find(node);
        auto it_after = neighbors_after.find(node);
        
        std::cout << "\nNode " << node << ":" << std::endl;
        std::cout << "  Position before: (" << points_before[node].x() << ", " 
                  << points_before[node].y() << ")" << std::endl;
        std::cout << "  Position after:  (" << points_after[node].x() << ", " 
                  << points_after[node].y() << ")" << std::endl;
        
        if (it_before != neighbors_before.end()) {
            std::cout << "  Neighbors before (" << it_before->second.size() << "): ";
            for (int nb : it_before->second) {
                std::cout << nb << " ";
            }
            std::cout << std::endl;
            
            std::cout << "    Distances: ";
            for (int nb : it_before->second) {
                double dx = points_before[nb].x() - points_before[node].x();
                double dy = points_before[nb].y() - points_before[node].y();
                double dist = std::sqrt(dx*dx + dy*dy);
                std::cout << dist << " ";
            }
            std::cout << std::endl;
        } else {
            std::cout << "  No neighbors before!" << std::endl;
        }
        
        if (it_after != neighbors_after.end()) {
            std::cout << "  Neighbors after  (" << it_after->second.size() << "): ";
            for (int nb : it_after->second) {
                std::cout << nb << " ";
            }
            std::cout << std::endl;
            
            std::cout << "    Distances: ";
            for (int nb : it_after->second) {
                double dx = points_after[nb].x() - points_after[node].x();
                double dy = points_after[nb].y() - points_after[node].y();
                double dist = std::sqrt(dx*dx + dy*dy);
                std::cout << dist << " ";
            }
            std::cout << std::endl;
        } else {
            std::cout << "  No neighbors after!" << std::endl;
        }
        
        if (it_before != neighbors_before.end() && it_after != neighbors_after.end()) {
            if (it_before->second == it_after->second) {
                std::cout << "  ✓ Neighbors UNCHANGED" << std::endl;
            } else {
                std::cout << "  ✗ Neighbors CHANGED" << std::endl;
            }
        }
    }
    std::cout << "\n========================\n" << std::endl;
}

void diagnosePointDistribution(const std::vector<Point2D>& points) {
    std::cout << "\n=== POINT DISTRIBUTION ANALYSIS ===" << std::endl;
    std::cout << "Total points: " << points.size() << std::endl;
    
    if (points.size() <= 20) {
        std::cout << "\nAll point positions:" << std::endl;
        for (size_t i = 0; i < points.size(); ++i) {
            std::cout << "  Node " << i << ": (" << points[i].x() << ", " 
                      << points[i].y() << ")" << std::endl;
        }
    }
    
    if (points.size() > 1) {
        std::vector<std::pair<double, int>> distances;
        for (size_t i = 0; i < points.size(); ++i) {
            if (i == 1) continue;
            double dx = points[i].x() - points[1].x();
            double dy = points[i].y() - points[1].y();
            double dist = std::sqrt(dx*dx + dy*dy);
            distances.push_back({dist, i});
        }
        
        std::sort(distances.begin(), distances.end());
        
        std::cout << "\nDistances from node 1 to all other nodes (sorted):" << std::endl;
        for (size_t i = 0; i < std::min(size_t(30), distances.size()); ++i) {
            std::cout << "  " << (i+1) << ". Node " << distances[i].second 
                      << ": distance = " << distances[i].first << std::endl;
        }
    }
    
    if (points.size() > 0) {
        double min_x = points[0].x(), max_x = points[0].x();
        double min_y = points[0].y(), max_y = points[0].y();
        
        for (const auto& p : points) {
            min_x = std::min(min_x, p.x());
            max_x = std::max(max_x, p.x());
            min_y = std::min(min_y, p.y());
            max_y = std::max(max_y, p.y());
        }
        
        std::cout << "\nBounding box:" << std::endl;
        std::cout << "  X: [" << min_x << ", " << max_x << "]" << std::endl;
        std::cout << "  Y: [" << min_y << ", " << max_y << "]" << std::endl;
        std::cout << "  Width: " << (max_x - min_x) << std::endl;
        std::cout << "  Height: " << (max_y - min_y) << std::endl;
    }
    
    std::cout << "====================================\n" << std::endl;
}

NeighborChangeInfo NeighborAnalyzer::compareNeighborsWithTolerance(
    const std::vector<Point2D>& points_before,
    const std::vector<Point2D>& points_after,
    const std::map<int, std::set<int>>& neighbors_before,
    const std::map<int, std::set<int>>& neighbors_after,
    double distance_tolerance
) {
    NeighborChangeInfo info;
    
    std::set<int> all_nodes;
    for (const auto& [node, _] : neighbors_before) {
        all_nodes.insert(node);
    }
    for (const auto& [node, _] : neighbors_after) {
        all_nodes.insert(node);
    }
    
    for (int node : all_nodes) {
        auto it_before = neighbors_before.find(node);
        auto it_after = neighbors_after.find(node);
        
        if (it_before == neighbors_before.end() || it_after == neighbors_after.end()) {
            continue;
        }
        
        const auto& before_set = it_before->second;
        const auto& after_set = it_after->second;
        
        // Si les ensembles sont identiques, pas de changement
        if (before_set == after_set) {
            continue;
        }
        
        // Vérifier si les voisins ajoutés/retirés sont à des distances similaires
        std::set<int> added, removed;
        std::set_difference(
            after_set.begin(), after_set.end(),
            before_set.begin(), before_set.end(),
            std::inserter(added, added.begin())
        );
        std::set_difference(
            before_set.begin(), before_set.end(),
            after_set.begin(), after_set.end(),
            std::inserter(removed, removed.begin())
        );
        
        // Calculer les distances moyennes des voisins ajoutés et retirés
        double avg_dist_added = 0.0, avg_dist_removed = 0.0;
        
        for (int nb : added) {
            double dx = points_after[nb].x() - points_after[node].x();
            double dy = points_after[nb].y() - points_after[node].y();
            avg_dist_added += std::sqrt(dx*dx + dy*dy);
        }
        avg_dist_added /= added.size();
        
        for (int nb : removed) {
            double dx = points_before[nb].x() - points_before[node].x();
            double dy = points_before[nb].y() - points_before[node].y();
            avg_dist_removed += std::sqrt(dx*dx + dy*dy);
        }
        avg_dist_removed /= removed.size();
        
        // Si les distances moyennes sont similaires, c'est probablement un échange de voisins équidistants
        double relative_diff = std::abs(avg_dist_added - avg_dist_removed) / 
                               std::max(avg_dist_added, avg_dist_removed);
        
        if (relative_diff > distance_tolerance) {
            // Vrai changement topologique
            info.has_changed = true;
            info.total_nodes_changed++;
            info.total_connections_added += added.size();
            info.total_connections_removed += removed.size();
            info.added_neighbors[node] = added;
            info.removed_neighbors[node] = removed;
            info.neighbor_count_change[node] = static_cast<int>(after_set.size()) - 
                                                static_cast<int>(before_set.size());
        }
    }
    
    return info;
}