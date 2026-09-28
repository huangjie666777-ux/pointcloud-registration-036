#include "icp3d/kdtree.h"

#include <algorithm>
#include <numeric>

namespace icp3d {

KdTree::KdTree(const PointCloud &points) : points_(points) {
  if (points_.empty()) {
    return;
  }
  std::vector<std::size_t> order(points_.size());
  std::iota(order.begin(), order.end(), std::size_t{0});
  nodes_.reserve(points_.size());
  build(order, 0, order.size(), 0);
}

int KdTree::build(std::vector<std::size_t> &order, std::size_t begin,
                  std::size_t end, int depth) {
  if (begin >= end) {
    return -1;
  }
  const int axis = depth % 3;
  const std::size_t mid = begin + (end - begin) / 2;
  auto first = order.begin() + static_cast<std::ptrdiff_t>(begin);
  auto middle = order.begin() + static_cast<std::ptrdiff_t>(mid);
  auto last = order.begin() + static_cast<std::ptrdiff_t>(end);
  std::nth_element(first, middle, last, [&](std::size_t a, std::size_t b) {
    return points_[a][axis] < points_[b][axis];
  });

  Node node;
  node.pointIndex = order[mid];
  const int self = static_cast<int>(nodes_.size());
  nodes_.push_back(node);
  nodes_[self].left = build(order, begin, mid, depth + 1);
  nodes_[self].right = build(order, mid + 1, end, depth + 1);
  return self;
}

KdTree::NearestResult KdTree::nearest(const Point &query) const {
  NearestResult best;
  if (nodes_.empty()) {
    return best;
  }
  best.index = nodes_[0].pointIndex;
  best.squaredDistance = (points_[best.index] - query).squaredNorm();
  search(0, query, 0, best);
  return best;
}

void KdTree::search(int nodeIndex, const Point &query, int depth,
                    NearestResult &best) const {
  if (nodeIndex < 0) {
    return;
  }
  const Node &node = nodes_[static_cast<std::size_t>(nodeIndex)];
  const Point &candidate = points_[node.pointIndex];
  const double dist2 = (candidate - query).squaredNorm();
  // 等距时取目标原始索引较小者，保证确定性。
  if (dist2 < best.squaredDistance ||
      (dist2 == best.squaredDistance && node.pointIndex < best.index)) {
    best.squaredDistance = dist2;
    best.index = node.pointIndex;
  }

  const int axis = depth % 3;
  const double delta = query[axis] - candidate[axis];
  const int nearSide = delta < 0.0 ? node.left : node.right;
  const int farSide = delta < 0.0 ? node.right : node.left;

  search(nearSide, query, depth + 1, best);
  if (delta * delta <= best.squaredDistance) {
    search(farSide, query, depth + 1, best);
  }
}

} // namespace icp3d
