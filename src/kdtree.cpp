#include "icp_registration/kdtree.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace icp {
namespace {

double coordinate(const Point& p, int axis) {
    return p[axis];
}

}  // namespace

KDTree::KDTree(const PointCloud& points) : points_(points) {
    order_.resize(points.size());
    for (std::size_t i = 0; i < points.size(); ++i) {
        order_[i] = i;
    }
    nodes_.reserve(points.size());
    if (!points.empty()) {
        build(0, order_.size(), 0);
    }
}

int KDTree::build(std::size_t begin, std::size_t end, int depth) {
    std::size_t count = end - begin;
    const int axis = depth % 3;

    // Split on the axis of largest spread of this subset for balance quality.
    double min_v = std::numeric_limits<double>::infinity();
    double max_v = -std::numeric_limits<double>::infinity();
    for (std::size_t i = begin; i < end; ++i) {
        const double v = coordinate(points_[order_[i]], axis);
        min_v = std::min(min_v, v);
        max_v = std::max(max_v, v);
    }
    double span_axis = max_v - min_v;
    int chosen_axis = axis;
    for (int candidate = 0; candidate < 3; ++candidate) {
        double lo = std::numeric_limits<double>::infinity();
        double hi = -std::numeric_limits<double>::infinity();
        for (std::size_t i = begin; i < end; ++i) {
            const double v = coordinate(points_[order_[i]], candidate);
            lo = std::min(lo, v);
            hi = std::max(hi, v);
        }
        if (hi - lo > span_axis) {
            span_axis = hi - lo;
            chosen_axis = candidate;
        }
    }

    std::sort(order_.begin() + static_cast<std::ptrdiff_t>(begin),
              order_.begin() + static_cast<std::ptrdiff_t>(end),
              [&](std::size_t a, std::size_t b) {
                  const double va = coordinate(points_[a], chosen_axis);
                  const double vb = coordinate(points_[b], chosen_axis);
                  if (va != vb) return va < vb;
                  return a < b;
              });

    const std::size_t median = begin + count / 2;
    Node node;
    node.index = order_[median];
    node.split_axis = chosen_axis;
    node.left = (median > begin) ? build(begin, median, depth + 1) : -1;
    node.right = (median + 1 < end)
                     ? build(median + 1, end, depth + 1)
                     : -1;
    nodes_.push_back(std::move(node));
    return static_cast<int>(nodes_.size() - 1);
}

std::size_t KDTree::nearestNeighbor(const Point& query,
                                    double& best_sq_dist) const {
    int best_node = -1;
    best_sq_dist = std::numeric_limits<double>::infinity();
    if (!nodes_.empty()) {
        // Root is the last node inserted (post-order build).
        search(static_cast<int>(nodes_.size() - 1), query, best_node,
               best_sq_dist);
    }
    return nodes_[static_cast<std::size_t>(best_node)].index;
}

void KDTree::search(int node_id,
                    const Point& query,
                    int& best_node,
                    double& best_sq_dist) const {
    if (node_id < 0) return;
    const Node& node = nodes_[static_cast<std::size_t>(node_id)];
    const Point& candidate = points_[node.index];
    const double dx = candidate.x() - query.x();
    const double dy = candidate.y() - query.y();
    const double dz = candidate.z() - query.z();
    const double sq = dx * dx + dy * dy + dz * dz;

    // Exact ties resolve to the smaller original target index.
    if (sq < best_sq_dist ||
        (sq == best_sq_dist &&
         (best_node < 0 || node.index < nodes_[static_cast<std::size_t>(best_node)].index))) {
        best_sq_dist = sq;
        best_node = node_id;
    }

    const double delta = query[node.split_axis] - candidate[node.split_axis];
    const double plane_sq = delta * delta;
    const int near_child = delta <= 0.0 ? node.left : node.right;
    const int far_child = delta <= 0.0 ? node.right : node.left;

    search(near_child, query, best_node, best_sq_dist);
    // Non-strict pruning: an equal-distance tie may live on the far side.
    if (plane_sq <= best_sq_dist) {
        search(far_child, query, best_node, best_sq_dist);
    }
}

}  // namespace icp
