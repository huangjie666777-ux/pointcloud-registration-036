#pragma once

#include "icp_registration/types.h"

namespace icp {

struct RigidEstimate {
    Transform transform = Transform::Identity();
    bool valid = false;
    bool degenerate = false;  // enough points, but geometry cannot fix R
};

// Throws std::invalid_argument on null, empty or non-finite clouds.
void validateFiniteCloud(const PointCloud& cloud, const char* name);

// Finite entries, R orthogonal and det(R)=+1, within fixed tolerances.
// Throws std::invalid_argument otherwise.
void validateRigidTransform(const Transform& transform, const char* name);

void validateOptions(const IcpOptions& options);

// Least-squares rigid transform mapping source[i] to target[i] (SVD,
// Arun/Horn/Kabsch/Umeyama). No scale, no reflection. Needs at least
// min_correspondences (>=3) pairs. Degenerate geometry (collinear/coincident
// centroids configuration that cannot determine the rotation, or a
// reflection-only alignment) is reported, never replaced by a guess.
RigidEstimate estimateRigidTransform(const PointCloud& source,
                                     const PointCloud& target,
                                     int min_correspondences = 3);

// Rotation angle in radians via the clamped acos trace formula.
double rotationAngle(const RotationMatrix& rotation);

}  // namespace icp
