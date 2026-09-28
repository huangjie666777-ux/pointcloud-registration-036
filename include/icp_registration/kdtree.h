#pragma once

#include "icp_registration/types.h"

#include <Eigen/Core>
#include <cstddef>

namespace icp {

// Static balanced 3D KD-tree over a fixed target cloud.
// The cloud must outlive the tree. No distance matrix is stored: memory is
// O(n) for the point-index permutation, and each query is O(log n) average.
class KDTree {
public:
    explicit KDTree(const PointCloud& points);

    // Nearest-neighbour query. On exact distance ties the neighbour with the
    // smaller original point index is returned.
    // Returns the original index; best_sq_dist receives the squared distance.
    std::size_t nearestNeighbor(const Point& query, double& best_sq_dist) const;

    std::size_t size() const noexcept { return points_.size(); }

private:
    struct Node {
        std::size_t index = 0;  // original point index
        int split_axis = 0;
        int left = -1;
        int right = -1;
    };

    const PointCloud& points_;
    std::vector<Node> nodes_;
    std::vector<std::size_t> order_;

    int build(std::size_t begin, std::size_t end, int depth);
    void search(int node_id,
                const Point& query,
                int& best_node,
                double& best_sq_dist) const;
};

}  // namespace icp
