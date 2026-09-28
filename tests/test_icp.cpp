#include "icp_registration/icp.h"

#include <Eigen/Geometry>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace icp;

namespace {

int g_failures = 0;
int g_checks = 0;

void check(bool condition, const std::string& message) {
    ++g_checks;
    if (!condition) {
        ++g_failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}

PointCloud makeGrid(int nx, int ny, int nz, double spacing = 0.5) {
    PointCloud cloud;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                cloud.emplace_back(i * spacing, j * spacing, k * spacing);
            }
        }
    }
    return cloud;
}

void testKdTreeTieBreaksByOriginalIndex() {
    const PointCloud target = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {-1.0, 0.0, 0.0},  // index 2, equidistant from query as index 1
    };
    const KDTree tree(target);
    double sq = -1.0;
    const std::size_t index = tree.nearestNeighbor(Point(0.0, 0.0, 0.0), sq);
    check(index == 0, "nearest point index");
    check(std::abs(sq) < 1e-20, "nearest point distance");

    sq = -1.0;
    const std::size_t tie_index =
        tree.nearestNeighbor(Point(0.0, 1.0, 0.0), sq);
    // Indices 0 (sq=1), then 1 and 2 both at sq=2: not a tie of the best.
    check(tie_index == 0, "best point wins over equidistant runners-up");

    const PointCloud symmetric = {
        {1.0, 0.0, 0.0},
        {-1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, -1.0, 0.0},
    };
    const KDTree tree2(symmetric);
    sq = -1.0;
    const std::size_t symmetric_index =
        tree2.nearestNeighbor(Point(0.0, 0.0, 0.0), sq);
    // All four are exactly equidistant: smallest original index must win.
    check(symmetric_index == 0, "exact tie resolves to smallest index");
    check(std::abs(sq - 1.0) < 1e-20, "tie distance preserved");
}

void testRigidEstimateKnownMotion() {
    const PointCloud source = makeGrid(4, 4, 4, 0.4);
    const Eigen::AngleAxisd angle(0.35, Eigen::Vector3d(0.3, 0.5, 0.8).normalized());
    Transform truth = Transform::Identity();
    truth.linear() = angle.toRotationMatrix();
    truth.translation() = Point(0.7, -0.4, 0.25);

    PointCloud target;
    for (const Point& p : source) target.push_back(truth * p);

    const RigidEstimate estimate = estimateRigidTransform(source, target);
    check(estimate.valid && !estimate.degenerate, "known motion is solvable");
    check((estimate.transform.matrix() - truth.matrix()).cwiseAbs().maxCoeff() < 1e-10,
          "SVD recovers known rigid motion");
    const RotationMatrix r = estimate.transform.linear();
    check(std::abs(r.determinant() - 1.0) < 1e-12, "estimated det(R)=+1");
}

void testRigidEstimateRejectsDegenerate() {
    const PointCloud source = {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {3, 0, 0}};
    const PointCloud target = {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {3, 0, 0}};
    const RigidEstimate collinear = estimateRigidTransform(source, target);
    check(!collinear.valid && collinear.degenerate,
          "collinear pairs are degenerate");

    const PointCloud same_a = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
    const PointCloud same_b = {{2, 2, 2}, {2, 2, 2}, {2, 2, 2}};
    const RigidEstimate coincident = estimateRigidTransform(same_a, same_b);
    check(!coincident.valid && coincident.degenerate,
          "coincident points are degenerate");
}

void testValidation() {
    const PointCloud good = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    const PointCloud empty;
    bool threw = false;
    try { validateFiniteCloud(empty, "empty"); } catch (const std::invalid_argument&) { threw = true; }
    check(threw, "empty cloud rejected");

    PointCloud bad = good;
    bad[1].x() = std::numeric_limits<double>::infinity();
    threw = false;
    try { validateFiniteCloud(bad, "bad"); } catch (const std::invalid_argument&) { threw = true; }
    check(threw, "non-finite coordinate rejected");

    Transform mirrored = Transform::Identity();
    mirrored.linear()(0, 0) = -1.0;
    threw = false;
    try { validateRigidTransform(mirrored, "mirror"); }
    catch (const std::invalid_argument&) { threw = true; }
    check(threw, "det(R)=-1 rejected");

    IcpOptions options;
    options.max_iterations = 0;
    threw = false;
    try { validateOptions(options); } catch (const std::invalid_argument&) { threw = true; }
    check(threw, "max_iterations=0 rejected");
}

void testIcpRecoversKnownTransformWithOutliers() {
    PointCloud source = makeGrid(6, 6, 6, 0.3);
    const Eigen::AngleAxisd angle(0.22, Eigen::Vector3d(0.2, -0.7, 0.5).normalized());
    Transform truth = Transform::Identity();
    truth.linear() = angle.toRotationMatrix();
    truth.translation() = Point(0.3, -0.2, 0.15);

    PointCloud target;
    for (const Point& p : source) target.push_back(truth * p);
    // Outlier target points far from anything, must never attract matches.
    target.push_back(truth * source[10] + Point(50.0, -30.0, 20.0));
    target.push_back(truth * source[20] + Point(-40.0, 60.0, 10.0));

    // Initial guess: deliberately perturbed around the truth.
    const Eigen::AngleAxisd perturb(0.08, Eigen::Vector3d::UnitZ());
    Transform initial = Transform::Identity();
    initial.linear() = perturb.toRotationMatrix() * truth.linear();
    initial.translation() = truth.translation() + Point(0.05, 0.05, 0.0);

    const double before = [&] {
        const KDTree tree(target);
        double sum = 0.0;
        for (const Point& p : source) {
            double sq = 0.0;
            tree.nearestNeighbor(initial * p, sq);
            sum += sq;
        }
        return std::sqrt(sum / static_cast<double>(source.size()));
    }();

    IcpOptions options;
    options.max_iterations = 100;
    options.max_correspondence_distance = 0.5;
    options.translation_tolerance = 1e-8;
    options.rotation_tolerance_rad = 1e-8;
    const IcpResult result = runIcp(source, target, initial, options);

    check(result.converged, "ICP converges for the known transform");
    check(result.reason == TerminationReason::Converged, "reason is converged");
    check(result.num_correspondences == source.size(),
          "every source point matched (outliers ignored)");
    check(result.rms_distance < 1e-7, "final RMS is near zero");
    check(before > result.rms_distance, "RMS improved over the initial guess");
    check((result.transform.matrix() - truth.matrix()).cwiseAbs().maxCoeff() < 1e-6,
          "final pose matches ground truth");
    check(std::abs(result.transform.linear().determinant() - 1.0) < 1e-12,
          "accumulated rotation stays proper");
}

void testIcpMaxIterationsReturnsUnconverged() {
    const PointCloud source = makeGrid(5, 5, 5, 0.4);
    const Eigen::AngleAxisd angle(0.5, Eigen::Vector3d::UnitY());
    Transform truth = Transform::Identity();
    truth.linear() = angle.toRotationMatrix();
    truth.translation() = Point(0.6, 0.1, -0.3);
    PointCloud target;
    for (const Point& p : source) target.push_back(truth * p);

    IcpOptions options;
    options.max_iterations = 2;
    options.max_correspondence_distance = 5.0;
    options.translation_tolerance = 1e-12;
    options.rotation_tolerance_rad = 1e-12;
    const IcpResult result =
        runIcp(source, target, Transform::Identity(), options);
    check(!result.converged, "tight thresholds with 2 iterations do not converge");
    check(result.reason == TerminationReason::MaxIterationsReached,
          "max-iterations termination reported");
    check(result.iterations == 2, "two iterations consumed");
    check(result.rms_distance >= 0.0 && std::isfinite(result.rms_distance),
          "final RMS still recomputed and returned");
}

void testIcpRejectsFarCorrespondences() {
    const PointCloud source = makeGrid(3, 3, 3, 0.5);
    PointCloud target = source;
    target.push_back(Point(100.0, 100.0, 100.0));  // extra unmatched target

    Transform shifted = Transform::Identity();
    shifted.translation() = Point(5.0, 0.0, 0.0);  // far from any target
    IcpOptions options;
    options.max_iterations = 10;
    options.max_correspondence_distance = 0.5;
    const IcpResult result =
        runIcp(source, target, shifted, options);
    check(!result.converged, "ICP cannot converge without correspondences");
    check(result.reason == TerminationReason::InsufficientCorrespondences,
          "insufficient correspondences reported");
    check(result.num_correspondences == 0, "no pair within the distance gate");
}

void testInputsAreNotModified() {
    const PointCloud source = makeGrid(3, 3, 3, 0.5);
    PointCloud target;
    Transform truth = Transform::Identity();
    truth.linear() = Eigen::AngleAxisd(0.1, Eigen::Vector3d::UnitX())
                         .toRotationMatrix();
    truth.translation() = Point(0.1, 0.0, 0.0);
    for (const Point& p : source) target.push_back(truth * p);

    const PointCloud source_copy = source;
    const PointCloud target_copy = target;
    Transform initial = Transform::Identity();
    IcpOptions options;
    runIcp(source, target, initial, options);
    check(source.size() == source_copy.size(), "source size unchanged");
    bool same = true;
    for (std::size_t i = 0; i < source.size(); ++i) {
        same &= (source[i] - source_copy[i]).norm() < 1e-15;
        same &= (target[i] - target_copy[i]).norm() < 1e-15;
    }
    check(same, "input cloud coordinates unchanged");
}

}  // namespace

int main() {
    testKdTreeTieBreaksByOriginalIndex();
    testRigidEstimateKnownMotion();
    testRigidEstimateRejectsDegenerate();
    testValidation();
    testIcpRecoversKnownTransformWithOutliers();
    testIcpMaxIterationsReturnsUnconverged();
    testIcpRejectsFarCorrespondences();
    testInputsAreNotModified();

    std::cout << g_checks - g_failures << '/' << g_checks << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
