#include <Eigen/Dense>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>

#include "icp3d/icp.h"
#include "icp3d/kdtree.h"

using namespace icp3d;

namespace {

double nearestRms(const PointCloud &source, const PointCloud &target,
                  const Eigen::Matrix4d &transform) {
  KdTree tree(target);
  double sum2 = 0.0;
  for (const auto &point : source) {
    Eigen::Vector4d homogeneous(point.x(), point.y(), point.z(), 1.0);
    const Eigen::Vector3d transformed =
        (transform * homogeneous).head<3>();
    sum2 += tree.nearest(transformed).squaredDistance;
  }
  return std::sqrt(sum2 / static_cast<double>(source.size()));
}

void printTransform(const Eigen::Matrix4d &transform) {
  Eigen::IOFormat format(7, 0, " ", "\n", "        [", "]", "[", "]");
  std::cout << transform.format(format) << "\n";
}

} // namespace

int main() {
  std::cout << std::setprecision(7);

  // 1) 生成目标点云：非平面、非退化的半椭球表面采样。
  PointCloud target;
  std::mt19937 rng(20260928);
  std::uniform_real_distribution<double> angle(0.0, 2.0 * M_PI);
  std::uniform_real_distribution<double> zcoord(-1.0, 1.0);
  for (int i = 0; i < 240; ++i) {
    const double theta = angle(rng);
    const double z = zcoord(rng);
    const double ring = std::sqrt(std::max(0.0, 1.0 - z * z));
    target.emplace_back(2.0 * ring * std::cos(theta) + 1.0,
                        ring * std::sin(theta) - 0.5,
                        1.5 * z + 0.3);
  }

  // 2) 已知真实“源->目标”刚体变换；源点取目标的子集并反变换。
  const double angleRad = 0.35;
  Eigen::Matrix3d trueRotation;
  trueRotation =
      Eigen::AngleAxisd(angleRad, Eigen::Vector3d(0.3, -0.8, 0.5).normalized());
  const Eigen::Vector3d trueTranslation(0.6, -0.4, 0.25);
  Eigen::Matrix4d trueTransform = Eigen::Matrix4d::Identity();
  trueTransform.topLeftCorner<3, 3>() = trueRotation;
  trueTransform.topRightCorner<3, 1>() = trueTranslation;
  const Eigen::Matrix4d invTrue = trueTransform.inverse();

  PointCloud source;
  for (std::size_t i = 0; i < target.size(); i += 3) {
    Eigen::Vector4d homogeneous(target[i].x(), target[i].y(),
                                target[i].z(), 1.0);
    source.push_back((invTrue * homogeneous).head<3>());
  }

  // 3) 加入离群源点：真实变换后远离任何目标点。
  source.emplace_back(50.0, -40.0, 30.0);
  source.emplace_back(-30.0, 40.0, -25.0);

  // 4) 初始姿态故意偏离真值（旋转约 0.14 rad，平移有误差）。
  Eigen::Matrix4d initial = Eigen::Matrix4d::Identity();
  initial.topLeftCorner<3, 3>() =
      Eigen::AngleAxisd(0.14, Eigen::Vector3d::UnitZ()).toRotationMatrix();
  initial.topRightCorner<3, 1>() = Eigen::Vector3d(0.25, -0.1, 0.1);

  IcpConfig config;
  config.maxIterations = 100;
  config.maxCorrespondenceDistance = 0.5;
  config.translationTolerance = 1e-8;
  config.rotationTolerance = 1e-8;

  std::cout << "=== 工件扫描点到点 ICP 配准示例 ===\n";
  std::cout << "源点数(含2个离群点): " << source.size()
            << ", 目标点数: " << target.size() << "\n";
  std::cout << "配准前最近点 RMS: " << nearestRms(source, target, initial)
            << "\n\n";

  const IcpResult result = alignPointToPoint(source, target, initial, config);

  std::cout << "终止原因: " << toString(result.reason)
            << (result.converged ? "（已收敛）" : "（未收敛）") << "\n";
  std::cout << "迭代次数: " << result.iterations
            << ", 最终有效配对数: " << result.numCorrespondences
            << "（离群点已被距离阈值剔除）\n";
  std::cout << "配准后最近点 RMS（最终变换重新匹配）: " << result.rms
            << "\n\n";
  std::cout << "最终源->目标变换:\n";
  printTransform(result.transform);
  std::cout << "\n真实变换（参考）:\n";
  printTransform(trueTransform);
  std::cout << "\n与真实变换的旋转角误差(rad): "
            << Eigen::AngleAxisd(
                   result.transform.topLeftCorner<3, 3>().transpose() *
                       trueTransform.topLeftCorner<3, 3>())
                   .angle()
            << "\n平移误差范数: "
            << (result.transform.topRightCorner<3, 1>() - trueTranslation)
                   .norm()
            << "\n";

  return result.converged ? 0 : 1;
}
