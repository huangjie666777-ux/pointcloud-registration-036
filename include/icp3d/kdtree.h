#pragma once

#include <Eigen/Dense>

#include <cstddef>
#include <limits>

#include "icp3d/types.h"

namespace icp3d {

// 目标点云的静态平衡 KD 树（3 维，双精度）。
// 只构建一次并在各轮 ICP 中复用；不保存点对距离矩阵。
class KdTree {
public:
  explicit KdTree(const PointCloud &points);

  struct NearestResult {
    std::size_t index = 0;
    double squaredDistance = std::numeric_limits<double>::infinity();
  };

  // 返回最近点在目标点云中的原始索引及其平方距离。
  // 多个点等距时，返回原始索引最小者。points 非空时结果有效。
  NearestResult nearest(const Point &query) const;

  std::size_t size() const { return points_.size(); }
  const PointCloud &points() const { return points_; }

private:
  struct Node {
    std::size_t pointIndex = 0;
    int left = -1;
    int right = -1;
  };

  const PointCloud &points_;
  std::vector<Node> nodes_;

  int build(std::vector<std::size_t> &order, std::size_t begin,
            std::size_t end, int depth);
  void search(int nodeIndex, const Point &query, int depth,
              NearestResult &best) const;
};

} // namespace icp3d
