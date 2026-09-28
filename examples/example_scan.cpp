#include "icp_registration/icp.h"

#include <Eigen/Geometry>
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace icp;

namespace {

// Synthetic workpiece surface: jittered points sampled on a solid box patch.
PointCloud makeWorkpieceScan(int nx, int ny, int nz, double spacing,
                             double jitter, unsigned int seed) {
    PointCloud cloud;
    unsigned int state = seed;
    auto next_noise = [&]() {
        state = state * 1664525u + 1013904223u;
        return (static_cast<double>(state % 2001) / 1000.0 - 1.0) * jitter;
    };
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int layer = 0; layer < nz; ++layer) {
                cloud.emplace_back(i * spacing + next_noise(),
                                   j * spacing + next_noise(),
                                   layer * spacing + next_noise());
            }
        }
    }
    return cloud;
}

void printPose(const char* label, const Transform& transform) {
    std::cout << label << "\n";
    const Eigen::IOFormat fmt(7, 0, " ", "\n", "    ", "", "", "");
    std::cout << std::fixed << std::setprecision(7);
    std::cout << "  R =\n    "
              << transform.linear().format(fmt) << "\n";
    std::cout << "  t = (" << transform.translation().x() << ", "
              << transform.translation().y() << ", "
              << transform.translation().z() << ")\n";
}

double rmsUnder(const PointCloud& source, const PointCloud& target,
                const Transform& transform) {
    const KDTree tree(target);
    double sum_sq = 0.0;
    for (const Point& point : source) {
        double sq = 0.0;
        tree.nearestNeighbor(transform * point, sq);
        sum_sq += sq;
    }
    return std::sqrt(sum_sq / static_cast<double>(source.size()));
}

}  // namespace

int main() {
    // Ground-truth placement of the second scan relative to the first.
    const Eigen::AngleAxisd truth_angle(
        0.18, Eigen::Vector3d(0.35, -0.55, 0.75).normalized());
    Transform truth = Transform::Identity();
    truth.linear() = truth_angle.toRotationMatrix();
    truth.translation() = Point(0.45, -0.30, 0.20);

    PointCloud scan_source = makeWorkpieceScan(8, 8, 8, 0.10, 0.001, 12345u);
    PointCloud scan_target;
    scan_target.reserve(scan_source.size() + 25);
    // The second scan has its own independent measurement noise.
    unsigned int target_state = 67890u;
    auto target_noise = [&]() {
        target_state = target_state * 1664525u + 1013904223u;
        return (static_cast<double>(target_state % 2001) / 1000.0 - 1.0) * 0.001;
    };
    for (const Point& point : scan_source) {
        Point measured = truth * point;
        measured += Point(target_noise(), target_noise(), target_noise());
        scan_target.push_back(measured);
    }
    // Non-overlapping outlier returns (clutter around the workpiece).
    scan_target.emplace_back(8.0, 6.0, 5.0);
    scan_target.emplace_back(-7.0, 9.0, -4.0);
    scan_target.emplace_back(10.0, -8.0, 6.0);

    // Coarse initial guess from scanner odometry: slightly off the truth.
    Transform initial = Transform::Identity();
    const Eigen::AngleAxisd perturb(0.04, Eigen::Vector3d::UnitZ());
    initial.linear() = perturb.toRotationMatrix() * truth.linear();
    initial.translation() = truth.translation() + Point(0.03, 0.02, -0.02);

    IcpOptions options;
    options.max_iterations = 80;
    options.max_correspondence_distance = 0.10;
    options.translation_tolerance = 1e-8;
    options.rotation_tolerance_rad = 1e-8;

    const double before_rms = rmsUnder(scan_source, scan_target, initial);
    const IcpResult result = runIcp(scan_source, scan_target, initial, options);

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "=== Workpiece scan registration example ===\n";
    std::cout << "source points : " << scan_source.size() << "\n";
    std::cout << "target points : " << scan_target.size()
              << " (incl. 3 clutter returns)\n\n";
    printPose("Ground truth (source -> target):", truth);
    printPose("\nInitial guess:", initial);
    printPose("\nICP result:", result.transform);

    const double pose_error =
        (result.transform.matrix() - truth.matrix()).cwiseAbs().maxCoeff();
    std::cout << "\n--- Metrics ---\n";
    std::cout << "RMS before alignment : " << before_rms << "\n";
    std::cout << "RMS after alignment  : " << result.rms_distance << "\n";
    std::cout << "iterations          : " << result.iterations << "\n";
    std::cout << "correspondences     : " << result.num_correspondences
              << " / " << scan_source.size() << "\n";
    std::cout << "pose error (max|.|) : " << pose_error << "\n";
    std::cout << "converged           : " << (result.converged ? "yes" : "no")
              << "\n";
    std::cout << "termination reason  : " << toString(result.reason) << "\n";
    return 0;
}
