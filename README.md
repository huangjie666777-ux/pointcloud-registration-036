# 工件扫描点云配准库 (`icp_registration`)

可复用的 C++17 三维点云刚体配准静态库，实现点到点 ICP（Iterative Closest Point），
用于把两次工件扫描对齐到同一坐标系。纯双精度、无缩放、无镜像；不承诺从任意初始位姿
收敛到全局最优（ICP 是局部方法，需要合理初始位姿与重叠区域）。

## 目录结构

- `include/icp_registration/types.h` — 公共数据类型、配准选项与结果、终止原因。
- `include/icp_registration/kdtree.h` / `src/kdtree.cpp` — 目标点云三维 KD 树空间索引。
- `include/icp_registration/rigid_transform.h` / `src/rigid_transform.cpp` —
  参数校验与 SVD 最小二乘刚体增量（Arun/Horn/Kabsch/Umeyama）。
- `include/icp_registration/icp.h` / `src/icp.cpp` — 对应点匹配与 ICP 迭代控制。
- `tests/test_icp.cpp` — 单元/集成测试（无第三方测试框架）。
- `examples/example_scan.cpp` — 已知变换 + 离群点的调用示例。
- `third_party/eigen3` — Eigen 3.4.0 头文件。

## 构建与运行

```sh
make            # 构建静态库、测试与示例
make test       # 运行测试
make example    # 运行含已知变换和离群点的示例
```

直接使用 g++：

```sh
g++ -std=c++17 -O2 -Iinclude -Ithird_party/eigen3 \
    src/*.cpp your_app.cpp -o your_app
```

## 快速使用

```cpp
#include "icp_registration/icp.h"

icp::PointCloud source = /* 扫描1，单位一致的双精度坐标 */;
icp::PointCloud target = /* 扫描2 */;
icp::Transform initial = icp::Transform::Identity(); /* 源->目标初始位姿 */

icp::IcpOptions options;
options.max_iterations = 80;
options.max_correspondence_distance = 0.15;   // 对应点最大距离（与坐标同单位）
options.translation_tolerance = 1e-8;        // 平移收敛阈值（单位同坐标）
options.rotation_tolerance_rad = 1e-8;       // 旋转收敛阈值（弧度）

icp::IcpResult result = icp::runIcp(source, target, initial, options);
// result.transform / iterations / num_correspondences / rms_distance
// result.converged / result.reason
```

## 算法说明

1. 目标 KD 树只构建一次并在所有迭代复用；不保存完整距离矩阵，存储 O(n)，
   单次最近邻平均 O(log n)。多个源点可匹配同一目标点；距离完全相等时，选择
   目标点原始索引较小者（搜索时不做等距剪枝以保证该规则）。
2. 每轮先用当前累计变换变换全部源点，再在目标中找最近点；平方距离超过
   `max_correspondence_distance^2` 的配对被剔除。
3. 对保留配对求质心，中心化构造 3x3 交叉协方差并做 JacobiSVD：
   `R = V U^T`；`det(R) < 0` 时按 Umeyama 方式在零奇异方向修正，禁止镜像。
   增量（当前帧）左乘到累计变换：`T_new = dT * T_old`，并对累计旋转做一次
   `U V^T` 重正交化，抑制舍入漂移，行列式始终为 +1。
4. 平移增量 `|t_new - t_old|` 与旋转增量 `angle(R_new R_old^T)`（acos 迹公式）
   同时不超过阈值才判定收敛；到达 `max_iterations` 仍未达标则返回
   `MaxIterationsReached`、`converged=false` 与当前估计。
5. 返回的配对数与 RMS 始终用**最终变换重新匹配**目标计算，不复用上一轮配对。

## 失败与校验语义

- 空点云、非有限（NaN/Inf）坐标抛 `std::invalid_argument`；输入不被修改。
- 初始变换要求有限、齐次末行 `[0 0 0 1]`、`R^T R = I`（偏差 ≤ 1e-10）且
  `det(R) = +1`（误差 ≤ 1e-10）；缩放/镜像一律拒绝。
- 选项要求 `max_iterations >= 1`、距离阈值为正有限值、收敛阈值非负有限、
  `min_correspondences >= 3`。
- 有效配对少于 `min_correspondences`（默认 3）返回
  `InsufficientCorrespondences`。
- SVD 中交叉矩阵相对最大奇异值的第二大奇异值 ≤ 1e-12（即共线/重合导致
  旋转不可确定），或满秩数据只能由镜像对齐时，返回 `DegenerateGeometry`，
  不会任取一个变换冒充成功。平面（秩 2）配置可以确定旋转，属合法情形。
- 数值容差汇总：正交性/行列式 1e-10；SVD 退化判定 1e-12（相对尺度）；
  收敛阈值由调用方给定，示例取 1e-8。
