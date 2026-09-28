#include "icp3d/icp.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "icp3d/kdtree.h"
#include "icp3d/rigid.h"

namespace icp3d {

namespace {

void validateCloud(const PointCloud &cloud, const char *name) {
  if (cloud.empty()) {
    throw std::invalid_argument(std::string(name) + "点云不能为空");
  }
  for (const Point &point : cloud) {
    if (!point.allFinite()) {
      throw std::invalid_argument(std::string(name) + "点云含非有限坐标");
    }
  }
}

void validateConfig(const IcpConfig &config) {
  if (config.maxIterations <= 0) {
    throw std::invalid_argument("maxIterations 必须为正整数");
  }
  if (!std::isfinite(config.maxCorrespondenceDistance) ||
      config.maxCorrespondenceDistance <= 0.0) {
    throw std::invalid_argument("maxCorrespondenceDistance 必须为正有限值");
  }
  if (!std::isfinite(config.translationTolerance) ||
      config.translationTolerance <= 0.0) {
    throw std::invalid_argument("translationTolerance 必须为正有限值");
  }
  if (!std::isfinite(config.rotationTolerance) ||
      config.rotationTolerance <= 0.0) {
    throw std::invalid_argument("rotationTolerance 必须为正有限值");
  }
  if (!std::isfinite(config.rankEpsilon) || config.rankEpsilon <= 0.0) {
    throw std::invalid_argument("rankEpsilon 必须为正有限值");
  }
}

Point applyTransform(const RigidTransform &pose, const Point &point) {
  return pose.rotation * point + pose.translation;
}

// 用当前姿态变换源点，在目标树中匹配，按最大距离过滤。
std::vector<Correspondence> matchPairs(const PointCloud &source,
                                       const KdTree &targetTree,
                                       const RigidTransform &pose,
                                       double maxDistance) {
  std::vector<Correspondence> pairs;
  pairs.reserve(source.size());
  const double maxDist2 = maxDistance * maxDistance;
  const PointCloud &target = targetTree.points();
  for (const Point &sourcePoint : source) {
    const Point transformed = applyTransform(pose, sourcePoint);
    const KdTree::NearestResult nearest = targetTree.nearest(transformed);
    if (nearest.squaredDistance <= maxDist2) {
      // 注意：左侧为“当前变换后的源点”，因此求出的是增量 Δ。
      pairs.push_back({transformed, target[nearest.index]});
    }
  }
  return pairs;
}

double rmsOverPairs(const std::vector<Correspondence> &pairs) {
  if (pairs.empty()) {
    return std::numeric_limits<double>::infinity();
  }
  double sum2 = 0.0;
  for (const auto &pair : pairs) {
    sum2 += (pair.target - pair.source).squaredNorm();
  }
  return std::sqrt(sum2 / static_cast<double>(pairs.size()));
}

} // namespace

IcpResult alignPointToPoint(const PointCloud &source,
                            const PointCloud &target,
                            const Eigen::Matrix4d &initialTransform,
                            const IcpConfig &config) {
  validateCloud(source, "源");
  validateCloud(target, "目标");
  validateConfig(config);
  RigidTransform pose = fromMatrix4(initialTransform);

  // 目标空间索引只构建一次，各轮迭代复用。
  KdTree targetTree(target);

  IcpResult result;
  result.transform = initialTransform;

  for (int iteration = 1; iteration <= config.maxIterations; ++iteration) {
    const std::vector<Correspondence> pairs =
        matchPairs(source, targetTree, pose, config.maxCorrespondenceDistance);
    if (pairs.size() < 3) {
      result.iterations = iteration;
      result.numCorrespondences = pairs.size();
      result.rms = rmsOverPairs(pairs);
      result.converged = false;
      result.reason = TerminationReason::InsufficientCorrespondences;
      result.message = "有效配对不足 3 组";
      result.transform = toMatrix4(pose);
      return result;
    }

    RigidTransform increment;
    try {
      increment = estimateRigid(pairs, config.rankEpsilon);
    } catch (const std::runtime_error &error) {
      result.iterations = iteration;
      result.numCorrespondences = pairs.size();
      result.rms = rmsOverPairs(pairs);
      result.converged = false;
      result.reason = TerminationReason::DegenerateGeometry;
      result.message = error.what();
      result.transform = toMatrix4(pose);
      return result;
    }

    // 正确累加：T_new = Δ * T_old（均为源->目标）。
    pose = compose(increment, pose);

    const double deltaTranslation = increment.translation.norm();
    const double deltaRotation = rotationAngle(increment.rotation);

    if (deltaTranslation < config.translationTolerance &&
        deltaRotation < config.rotationTolerance) {
      // 指标按最终变换重新匹配计算，不混用上轮误差。
      const std::vector<Correspondence> finalPairs =
          matchPairs(source, targetTree, pose,
                     config.maxCorrespondenceDistance);
      result.iterations = iteration;
      result.numCorrespondences = finalPairs.size();
      result.rms = rmsOverPairs(finalPairs);
      result.converged = true;
      result.reason = TerminationReason::Converged;
      result.message = "平移与旋转增量均低于收敛阈值";
      result.transform = toMatrix4(pose);
      return result;
    }
  }

  // 达到次数上限：返回未收敛状态及当前估计，指标按当前估计重新匹配。
  const std::vector<Correspondence> finalPairs =
      matchPairs(source, targetTree, pose, config.maxCorrespondenceDistance);
  result.iterations = config.maxIterations;
  result.numCorrespondences = finalPairs.size();
  result.rms = rmsOverPairs(finalPairs);
  result.converged = false;
  result.reason = TerminationReason::MaxIterationsReached;
  result.message = "达到最大迭代次数仍未收敛";
  result.transform = toMatrix4(pose);
  return result;
}

} // namespace icp3d
