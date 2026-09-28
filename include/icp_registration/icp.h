#pragma once

#include "icp_registration/kdtree.h"
#include "icp_registration/rigid_transform.h"
#include "icp_registration/types.h"

namespace icp {

struct Correspondence {
    std::size_t source_index;
    std::size_t target_index;
    double squared_distance;
};

// Match every transformed source point to its nearest target point. Pairs
// farther than max_correspondence_distance are rejected. Multiple source
// points may share one target; exact ties resolve to the smaller original
// target index inside KDTree.
std::vector<Correspondence> findCorrespondences(const PointCloud& source_transformed,
                                                const KDTree& target_tree,
                                                double max_correspondence_distance);

// Point-to-point ICP. Inputs are never modified. The initial guess (and the
// returned transform) maps source coordinates into the target frame.
IcpResult runIcp(const PointCloud& source,
                 const PointCloud& target,
                 const Transform& initial_transform,
                 const IcpOptions& options);

}  // namespace icp
