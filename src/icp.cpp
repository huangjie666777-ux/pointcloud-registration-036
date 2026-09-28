#include "icp_registration/icp.h"

#include <cmath>

namespace icp {
namespace {

PointCloud transformCloud(const PointCloud& cloud, const Transform& transform) {
    PointCloud moved;
    moved.reserve(cloud.size());
    for (const Point& point : cloud) {
        moved.push_back(transform * point);
    }
    return moved;
}

}  // namespace

std::vector<Correspondence> findCorrespondences(
    const PointCloud& source_transformed,
    const KDTree& target_tree,
    double max_correspondence_distance) {
    std::vector<Correspondence> correspondences;
    correspondences.reserve(source_transformed.size());
    const double max_sq =
        max_correspondence_distance * max_correspondence_distance;

    for (std::size_t i = 0; i < source_transformed.size(); ++i) {
        double sq_dist = 0.0;
        const std::size_t target_index =
            target_tree.nearestNeighbor(source_transformed[i], sq_dist);
        if (sq_dist <= max_sq) {
            correspondences.push_back({i, target_index, sq_dist});
        }
    }
    return correspondences;
}

IcpResult runIcp(const PointCloud& source,
                 const PointCloud& target,
                 const Transform& initial_transform,
                 const IcpOptions& options) {
    validateFiniteCloud(source, "source");
    validateFiniteCloud(target, "target");
    validateRigidTransform(initial_transform, "initial_transform");
    validateOptions(options);

    IcpResult result;
    result.transform = initial_transform;

    // The target index is built once and reused every iteration.
    const KDTree target_tree(target);

    bool converged = false;
    int iteration = 0;
    for (; iteration < options.max_iterations; ++iteration) {
        // 1. Transform the source with the current accumulated estimate.
        const PointCloud moved = transformCloud(source, result.transform);

        // 2. Nearest-neighbour matching in the fixed target index.
        const std::vector<Correspondence> correspondences =
            findCorrespondences(moved, target_tree,
                                options.max_correspondence_distance);
        if (static_cast<int>(correspondences.size()) <
            options.min_correspondences) {
            result.iterations = iteration;
            result.reason = TerminationReason::InsufficientCorrespondences;
            return result;
        }

        // 3. Rigid least-squares increment for the retained pairs. The SVD
        //    solver operates in the *current* frame: moved point -> target.
        PointCloud paired_source;
        PointCloud paired_target;
        paired_source.reserve(correspondences.size());
        paired_target.reserve(correspondences.size());
        for (const Correspondence& pair : correspondences) {
            paired_source.push_back(moved[pair.source_index]);
            paired_target.push_back(target[pair.target_index]);
        }
        const RigidEstimate estimate = estimateRigidTransform(
            paired_source, paired_target, options.min_correspondences);
        if (!estimate.valid) {
            result.iterations = iteration;
            result.reason = estimate.degenerate
                                ? TerminationReason::DegenerateGeometry
                                : TerminationReason::InsufficientCorrespondences;
            return result;
        }

        // 4. Accumulate: total <- increment * total (both map source to target).
        const Transform previous = result.transform;
        result.transform = Transform(estimate.transform.matrix() *
                                     previous.matrix());
        // Renormalize against accumulated round-off; stays a proper rotation.
        Eigen::JacobiSVD<Eigen::Matrix3d> svd(
            result.transform.linear(),
            Eigen::ComputeFullU | Eigen::ComputeFullV);
        result.transform.linear() =
            svd.matrixU() * svd.matrixV().transpose();

        // 5. Increment size is measured on the pose update delta_total =
        //    increment = new * old^{-1}, comparing frames directly.
        const double delta_translation =
            (result.transform.translation() - previous.translation()).norm();
        const RotationMatrix delta_rotation_matrix =
            result.transform.linear() * previous.linear().transpose();
        const double delta_rotation = rotationAngle(delta_rotation_matrix);

        if (delta_translation <= options.translation_tolerance &&
            delta_rotation <= options.rotation_tolerance_rad) {
            converged = true;
            iteration += 1;
            break;
        }
    }

    result.iterations = iteration;

    // Final metrics: re-match with the final transform, never reuse old pairs.
    const PointCloud final_moved = transformCloud(source, result.transform);
    const std::vector<Correspondence> final_correspondences =
        findCorrespondences(final_moved, target_tree,
                            options.max_correspondence_distance);
    result.num_correspondences = final_correspondences.size();
    if (!final_correspondences.empty()) {
        double sum_sq = 0.0;
        for (const Correspondence& pair : final_correspondences) {
            sum_sq += pair.squared_distance;
        }
        result.rms_distance =
            std::sqrt(sum_sq / static_cast<double>(final_correspondences.size()));
    }

    result.converged = converged;
    result.reason = converged ? TerminationReason::Converged
                              : TerminationReason::MaxIterationsReached;
    return result;
}

}  // namespace icp
