#pragma once

#include <Eigen/Dense>

#include "icp3d/types.h"

namespace icp3d {

// 点到点 ICP：将 source 对齐到 target。
//
// 输入约定：
//  - source/target 非空，坐标全部有限，两点集大小可不同、无需预设对应；
//  - initialTransform 为源->目标初始刚体变换，旋转须正交且行列式为 1；
//  - 输入不会被修改。配置项为正时抛 std::invalid_argument。
//
// 每轮：先用当前变换变换源点，再在目标 KD 树上找最近点，多个源点可
// 匹配同一目标，等距时取目标原始索引较小者；距离超阈值的配对剔除。
//
// 有效配对不足 3 组或几何退化时结果 converged=false 并给出原因，
// 不会返回任意冒充的变换。达到次数上限返回未收敛及当前估计。
// 输出 RMS 用最终变换重新匹配计算，不沿用上一轮误差。
IcpResult alignPointToPoint(const PointCloud &source,
                            const PointCloud &target,
                            const Eigen::Matrix4d &initialTransform,
                            const IcpConfig &config);

} // namespace icp3d
