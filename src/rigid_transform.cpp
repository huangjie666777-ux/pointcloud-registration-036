#include "icp_registration/rigid_transform.h"

#include <Eigen/SVD>
#include <cmath>
#include <stdexcept>
#include <string>

namespace icp {
namespace {

constexpr double kOrthogonalityTol = 1e-10;
constexpr double kDetTol = 1e-10;
// Relative singular-value cutoff for affine rank / reflection decisions.
constexpr double kSingularValueTol = 1e-12;

bool isFinite(const Point& p) {
    return std::isfinite(p.x()) && std::isfinite(p.y()) &&
           std::isfinite(p.z());
}

}  // namespace

void validateFiniteCloud(const PointCloud& cloud, const char* name) {
    if (cloud.empty()) {
        throw std::invalid_argument(std::string(name) + " point cloud is empty");
    }
    for (std::size_t i = 0; i < cloud.size(); ++i) {
        if (!isFinite(cloud[i])) {
            throw std::invalid_argument(std::string(name) +
                                        " contains a non-finite coordinate at index " +
                                        std::to_string(i));
        }
    }
}

void validateRigidTransform(const Transform& transform, const char* name) {
    const auto matrix = transform.matrix();
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 4; ++col) {
            if (!std::isfinite(matrix(row, col))) {
                throw std::invalid_argument(std::string(name) +
                                            " contains a non-finite entry");
            }
        }
    }
    const Eigen::Vector4d expected_last(0.0, 0.0, 0.0, 1.0);
    if ((matrix.row(3).transpose() - expected_last).cwiseAbs().maxCoeff() >
        kOrthogonalityTol) {
        throw std::invalid_argument(std::string(name) +
                                    " is not a homogeneous rigid transform");
    }
    const RotationMatrix r = transform.linear();
    const RotationMatrix deviation = r.transpose() * r - RotationMatrix::Identity();
    if (deviation.cwiseAbs().maxCoeff() > kOrthogonalityTol) {
        throw std::invalid_argument(std::string(name) +
                                    " rotation is not orthogonal");
    }
    if (std::abs(r.determinant() - 1.0) > kDetTol) {
        throw std::invalid_argument(std::string(name) +
                                    " rotation determinant is not +1 (mirroring rejected)");
    }
}

void validateOptions(const IcpOptions& options) {
    if (options.max_iterations < 1) {
        throw std::invalid_argument("max_iterations must be >= 1");
    }
    if (!std::isfinite(options.max_correspondence_distance) ||
        options.max_correspondence_distance <= 0.0) {
        throw std::invalid_argument(
            "max_correspondence_distance must be finite and positive");
    }
    if (!std::isfinite(options.translation_tolerance) ||
        options.translation_tolerance < 0.0) {
        throw std::invalid_argument("translation_tolerance must be finite and >= 0");
    }
    if (!std::isfinite(options.rotation_tolerance_rad) ||
        options.rotation_tolerance_rad < 0.0) {
        throw std::invalid_argument(
            "rotation_tolerance_rad must be finite and >= 0");
    }
    if (options.min_correspondences < 3) {
        throw std::invalid_argument("min_correspondences must be >= 3");
    }
}

RigidEstimate estimateRigidTransform(const PointCloud& source,
                                     const PointCloud& target,
                                     int min_correspondences) {
    RigidEstimate result;
    const auto n = source.size();
    if (min_correspondences < 3 ||
        n != target.size() ||
        static_cast<int>(n) < min_correspondences) {
        return result;
    }

    Point source_centroid = Point::Zero();
    Point target_centroid = Point::Zero();
    for (std::size_t i = 0; i < n; ++i) {
        source_centroid += source[i];
        target_centroid += target[i];
    }
    source_centroid /= static_cast<double>(n);
    target_centroid /= static_cast<double>(n);

    Eigen::Matrix3d cross = Eigen::Matrix3d::Zero();
    for (std::size_t i = 0; i < n; ++i) {
        cross += (source[i] - source_centroid) *
                 (target[i] - target_centroid).transpose();
    }

    Eigen::JacobiSVD<Eigen::Matrix3d> svd(
        cross, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Eigen::Vector3d singular = svd.singularValues();

    // Relative scale of the centered data; the two largest singular values
    // must be nonzero for the rotation to be determined (planar rank 2 is
    // sufficient; rank 1/collinear or rank 0 are not).
    const double scale = std::max(singular(0), 1.0);
    if (singular(0) <= kSingularValueTol * scale ||
        singular(1) <= kSingularValueTol * scale) {
        result.degenerate = true;
        return result;
    }

    RotationMatrix rotation =
        svd.matrixV() * svd.matrixU().transpose();
    if (rotation.determinant() < 0.0) {
        // Proper rotation impossible without reflection: if the smallest
        // singular direction carries real information the data wants a
        // mirror, so report degeneracy instead of forcing an answer.
        if (singular(2) > kSingularValueTol * scale) {
            result.degenerate = true;
            return result;
        }
        Eigen::Matrix3d corrected_v = svd.matrixV();
        corrected_v.col(2) *= -1.0;
        rotation = corrected_v * svd.matrixU().transpose();
    }

    Transform transform = Transform::Identity();
    transform.linear() = rotation;
    transform.translation() = target_centroid - rotation * source_centroid;
    result.transform = transform;
    result.valid = true;
    return result;
}

double rotationAngle(const RotationMatrix& rotation) {
    const double trace = rotation.trace();
    double cos_angle = 0.5 * (trace - 1.0);
    cos_angle = std::max(-1.0, std::min(1.0, cos_angle));
    return std::acos(cos_angle);
}

const char* toString(TerminationReason reason) {
    switch (reason) {
        case TerminationReason::NotStarted:
            return "not started";
        case TerminationReason::Converged:
            return "converged";
        case TerminationReason::MaxIterationsReached:
            return "max iterations reached";
        case TerminationReason::InsufficientCorrespondences:
            return "insufficient correspondences";
        case TerminationReason::DegenerateGeometry:
            return "degenerate geometry";
        case TerminationReason::InvalidArgument:
            return "invalid argument";
    }
    return "unknown";
}

}  // namespace icp
