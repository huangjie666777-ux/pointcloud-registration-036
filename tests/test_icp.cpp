#include <Eigen/Dense>

#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "icp3d/icp.h"
#include "icp3d/kdtree.h"
#include "icp3d/rigid.h"

using namespace icp3d;

namespace {

int g_failures = 0;
int g_checks = 0;

void check(bool condition, const std::string &name) {
  ++g_checks;
  if (condition) {
    std::cout << "  [PASS] " << name << "\n";
  } else {
    ++g_failures;
    std::cout << "  [FAIL] " << name << "\n";
  }
}

void checkThrows(const std::function<void()> &fn, const std::string &name) {
  bool thrown = false;
  try {
    fn();
  } catch (const std::invalid_argument &) {
    thrown = true;
  } catch (const std::exception &) {
  }
  check(thrown, name);
}

Eigen::Matrix4d makeTransform(const Eigen::Matrix3d &rotation,
                              const Eigen::Vector3d &translation) {
  Eigen::Matrix4d transform = Eigen::Matrix4d::Identity();
  transform.topLeftCorner<3, 3>() = rotation;
  transform.topRightCorner<3, 1>() = translation;
  return transform;
}

void testKdTree() {
  std::cout << "[KdTree]\n";
  PointCloud points;
  points.emplace_back(0, 0, 0);   // 0
  points.emplace_back(2, 0, 0);   // 1
  points.emplace_back(-2, 0, 0);  // 2
  points.emplace_back(0, 3, 0);   // 3
  points.emplace_back(0, -3, 0);  // 4
  KdTree tree(points);
  const auto nearest = tree.nearest(Point(0.9, 0.0, 0.0));
  check(nearest.index == 0 && std::abs(nearest.squaredDistance - 0.81) < 1e-12,
        "最近点查询返回正确索引与平方距离");

  // 查询点 (1,0,0) 与索引 0、1 等距，应返回较小原始索引。
  const auto tie = tree.nearest(Point(1.0, 0.0, 0.0));
  check(tie.index == 0, "等距时返回目标原始索引较小者");
}

void testRigidEstimate() {
  std::cout << "[RigidEstimate]\n";
  const Eigen::Matrix3d rotation =
      Eigen::AngleAxisd(0.7, Eigen::Vector3d(1.0, 2.0, 3.0).normalized())
          .toRotationMatrix();
  const Eigen::Vector3d translation(1.5, -2.0, 0.4);

  PointCloud source{
      {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1}, {2, -1, 0.5}};
  std::vector<Correspondence> pairs;
  for (const auto &point : source) {
    pairs.push_back({point, rotation * point + translation});
  }
  const RigidTransform estimated = estimateRigid(pairs);
  check(estimated.rotation.isApprox(rotation, 1e-10), "SVD 恢复已知旋转");
  check(estimated.translation.isApprox(translation, 1e-10), "SVD 恢复已知平移");
  check(std::abs(estimated.rotation.determinant() - 1.0) < 1e-12,
        "估计旋转行列式为 1");
  check((estimated.rotation.transpose() * estimated.rotation -
         Eigen::Matrix3d::Identity())
            .norm() < 1e-12,
        "估计旋转保持正交");

  bool failedForFew = false;
  try {
    estimateRigid({pairs[0], pairs[1]});
  } catch (const std::runtime_error &) {
    failedForFew = true;
  }
  check(failedForFew, "配对不足 3 组时明确失败");

  // 共线点：旋转不可唯一确定，必须判为退化。
  std::vector<Correspondence> collinear;
  for (double t : {0.0, 1.0, 2.0, 3.0}) {
    Eigen::Vector3d point(t, 0.0, 0.0);
    collinear.push_back({point, rotation * point + translation});
  }
  bool failedForDegenerate = false;
  try {
    estimateRigid(collinear);
  } catch (const std::runtime_error &) {
    failedForDegenerate = true;
  }
  check(failedForDegenerate, "共线几何退化时报错而非任取变换");

  // compose: T_new = Δ * T_old。
  RigidTransform pose{RigidTransform{rotation, translation}};
  RigidTransform delta{RigidTransform{
      Eigen::AngleAxisd(0.2, Eigen::Vector3d::UnitY()).toRotationMatrix(),
      Eigen::Vector3d(0.1, 0.2, 0.3)}};
  const RigidTransform composed = compose(delta, pose);
  const auto m = toMatrix4(composed);
  const Eigen::Matrix4d expected = toMatrix4(delta) * toMatrix4(pose);
  check(m.isApprox(expected), "增量按 Δ*T 正确累加");
  check(toMatrix4(fromMatrix4(m)).isApprox(m), "齐次矩阵互逆转换一致");

  bool invalidRotation = false;
  try {
    Eigen::Matrix3d mirror = Eigen::Matrix3d::Identity();
    mirror(2, 2) = -1.0;
    validateRotation(mirror);
  } catch (const std::invalid_argument &) {
    invalidRotation = true;
  }
  check(invalidRotation, "镜像旋转被拒绝");

  bool scaledRotation = false;
  try {
    Eigen::Matrix3d scaled = rotation * 2.0;
    validateRotation(scaled);
  } catch (const std::invalid_argument &) {
    scaledRotation = true;
  }
  check(scaledRotation, "含缩放的旋转被拒绝");
}

PointCloud makeSphereCloud(std::size_t count) {
  std::mt19937 rng(7);
  std::uniform_real_distribution<double> uniform(-1.0, 1.0);
  PointCloud cloud;
  while (cloud.size() < count) {
    Point point(uniform(rng), uniform(rng), uniform(rng));
    if (point.norm() <= 1.0) {
      cloud.push_back(point * 1.5);
    }
  }
  return cloud;
}

void testIcp() {
  std::cout << "[ICP]\n";
  const PointCloud target = makeSphereCloud(300);
  const Eigen::Matrix3d rotation =
      Eigen::AngleAxisd(0.25, Eigen::Vector3d(1.0, -1.0, 1.0).normalized())
          .toRotationMatrix();
  const Eigen::Vector3d translation(0.3, -0.2, 0.15);
  const Eigen::Matrix4d truth = makeTransform(rotation, translation);
  const Eigen::Matrix4d truthInv = truth.inverse();

  PointCloud source;
  for (std::size_t i = 0; i < target.size(); i += 4) {
    Eigen::Vector4d homogeneous(target[i].x(), target[i].y(),
                                target[i].z(), 1.0);
    source.push_back((truthInv * homogeneous).head<3>());
  }
  const PointCloud sourceCopy = source;
  const PointCloud targetCopy = target;

  IcpConfig config;
  config.maxIterations = 100;
  config.maxCorrespondenceDistance = 0.5;
  config.translationTolerance = 1e-9;
  config.rotationTolerance = 1e-9;

  const Eigen::Matrix4d initial = Eigen::Matrix4d::Identity();
  const IcpResult result = alignPointToPoint(source, target, initial, config);
  check(result.converged, "已知变换场景下 ICP 收敛");
  check(result.reason == TerminationReason::Converged, "终止原因标记为 Converged");
  check(result.transform.isApprox(truth, 1e-6), "最终变换恢复已知刚体变换");
  check(result.rms < 1e-5, "配准后 RMS 接近零");
  check(result.numCorrespondences == source.size(),
        "所有正常源点形成有效配对");
  check(source == sourceCopy && target == targetCopy, "配准过程不修改输入点云");

  // 离群点：应被最大对应距离剔除，配准仍成功。
  PointCloud sourceWithOutliers = source;
  sourceWithOutliers.emplace_back(20.0, 20.0, 20.0);
  sourceWithOutliers.emplace_back(-20.0, -20.0, -20.0);
  const IcpResult resultOutliers =
      alignPointToPoint(sourceWithOutliers, target, initial, config);
  check(resultOutliers.converged &&
            resultOutliers.numCorrespondences == source.size(),
        "离群源点被剔除且不影响配准");

  // 阈值过小：无有效配对，明确失败而非冒充成功。
  IcpConfig noPairs = config;
  noPairs.maxCorrespondenceDistance = 1e-9;
  const IcpResult resultNoPairs =
      alignPointToPoint(source, target, initial, noPairs);
  check(!resultNoPairs.converged &&
            resultNoPairs.reason ==
                TerminationReason::InsufficientCorrespondences,
        "有效配对不足时返回 InsufficientCorrespondences");

  // 达到次数上限：未收敛，返回当前估计与重新匹配的指标。
  IcpConfig capped = config;
  capped.maxIterations = 1;
  capped.translationTolerance = 1e-15;
  capped.rotationTolerance = 1e-15;
  const IcpResult resultCapped =
      alignPointToPoint(source, target, initial, capped);
  check(!resultCapped.converged && resultCapped.iterations == 1 &&
            resultCapped.reason == TerminationReason::MaxIterationsReached,
        "达到次数上限返回未收敛及当前估计");
}

void testValidation() {
  std::cout << "[Validation]\n";
  PointCloud cloud;
  cloud.emplace_back(0, 0, 0);
  cloud.emplace_back(1, 1, 1);
  PointCloud bad = cloud;
  bad[0].x() = std::numeric_limits<double>::infinity();
  const auto identity = Eigen::Matrix4d::Identity();
  IcpConfig config;

  PointCloud empty;
  checkThrows([&] { alignPointToPoint(empty, cloud, identity, config); },
              "拒绝空源点云");
  checkThrows([&] { alignPointToPoint(cloud, empty, identity, config); },
              "拒绝空目标点云");
  checkThrows([&] { alignPointToPoint(bad, cloud, identity, config); },
              "拒绝非有限坐标");

  Eigen::Matrix4d mirror = Eigen::Matrix4d::Identity();
  mirror(2, 2) = -1.0;
  checkThrows([&] { alignPointToPoint(cloud, cloud, mirror, config); },
              "拒绝镜像初始变换");

  IcpConfig badIterations = config;
  badIterations.maxIterations = 0;
  checkThrows(
      [&] { alignPointToPoint(cloud, cloud, identity, badIterations); },
      "拒绝非正最大迭代次数");
  IcpConfig badDistance = config;
  badDistance.maxCorrespondenceDistance = -1.0;
  checkThrows(
      [&] { alignPointToPoint(cloud, cloud, identity, badDistance); },
      "拒绝非正最大对应距离");
  IcpConfig badTolerance = config;
  badTolerance.translationTolerance = 0.0;
  checkThrows(
      [&] { alignPointToPoint(cloud, cloud, identity, badTolerance); },
      "拒绝非正平移阈值");
}

} // namespace

int main() {
  testKdTree();
  testRigidEstimate();
  testIcp();
  testValidation();
  std::cout << "\n共 " << g_checks << " 项检查，失败 " << g_failures
            << " 项\n";
  return g_failures == 0 ? 0 : 1;
}
