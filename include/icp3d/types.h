#pragma once

#include <Eigen/Dense>

#include <cstddef>
#include <string>
#include <vector>

namespace icp3d {

using Point = Eigen::Vector3d;
using PointCloud = std::vector<Point>;

// ICP 迭代终止原因。
enum class TerminationReason {
  Converged,                 // 平移与旋转增量均小于阈值
  MaxIterationsReached,      // 达到最大迭代次数，仍未收敛
  InsufficientCorrespondences, // 有效配对不足 3 组
  DegenerateGeometry,        // 配对几何退化，无法唯一确定旋转
};

const char *toString(TerminationReason reason);

// 点到点 ICP 参数。所有参数必须为正的有限值（见 icp.h 校验规则）。
struct IcpConfig {
  int maxIterations = 50;            // 最大迭代轮数，必须 > 0
  double maxCorrespondenceDistance = 1.0; // 最大对应距离，超过则剔除
  double translationTolerance = 1e-6; // 平移增量收敛阈值（与点坐标同单位）
  double rotationTolerance = 1e-6;   // 旋转增量收敛阈值（弧度）
  // 旋转可解的奇异值判据：最小奇异值必须大于
  // rankEpsilon * max(1, 最大奇异值)。平面/共线点云会被判为退化。
  double rankEpsilon = 1e-10;
};

// 配准结果。transform 始终为最终（或失败时当前）源到目标刚体变换。
struct IcpResult {
  Eigen::Matrix4d transform = Eigen::Matrix4d::Identity();
  int iterations = 0;
  std::size_t numCorrespondences = 0;
  double rms = 0.0; // 用最终变换重新匹配、按距离阈值过滤后的均方根距离
  bool converged = false;
  TerminationReason reason = TerminationReason::MaxIterationsReached;
  std::string message;
};

} // namespace icp3d
