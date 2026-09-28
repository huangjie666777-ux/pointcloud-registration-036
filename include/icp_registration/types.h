#pragma once

#include <Eigen/Dense>
#include <cstddef>
#include <limits>
#include <vector>

namespace icp {

using Point = Eigen::Vector3d;
using PointCloud = std::vector<Point>;
using RotationMatrix = Eigen::Matrix3d;
using Transform = Eigen::Isometry3d;

enum class TerminationReason {
    NotStarted = 0,
    Converged = 1,
    MaxIterationsReached = 2,
    InsufficientCorrespondences = 3,
    DegenerateGeometry = 4,
    InvalidArgument = 5,
};

struct IcpOptions {
    int max_iterations = 50;
    double max_correspondence_distance = 1.0;
    double translation_tolerance = 1e-6;
    double rotation_tolerance_rad = 1e-6;
    int min_correspondences = 3;
};

struct IcpResult {
    Transform transform = Transform::Identity();
    int iterations = 0;
    std::size_t num_correspondences = 0;
    double rms_distance = std::numeric_limits<double>::quiet_NaN();
    bool converged = false;
    TerminationReason reason = TerminationReason::NotStarted;
};

const char* toString(TerminationReason reason);

}  // namespace icp
