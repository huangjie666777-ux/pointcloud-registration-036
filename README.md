# 工件扫描点云刚体配准库 `icp3d`

可复用的 C++17 三维点云**点到点 ICP（Iterative Closest Point）**刚体配准库，
用于将两次工件扫描对齐到同一坐标系。基于 Eigen 3.4.0，仅含头文件依赖，
构建为静态库 `libicp3d.a`。

## 功能范围

- 输入源点云、目标点云与源到目标的初始刚体变换；双精度、同单位坐标。
- 两点集大小可不同，无需预设对应关系；输入不会被修改。
- 只求解刚体变换（旋转 + 平移），**不求解缩放、不接受镜像**。
- ICP 是局部方法：**不承诺从任意初始位姿收敛到全局最优**，需提供合理初值。

## 目录结构

```
include/icp3d/
  types.h    公共类型：Point/PointCloud、IcpConfig、IcpResult、终止原因
  kdtree.h   目标点云 3D KD 树（空间检索职责）
  rigid.h    刚体变换校验、组合与最小二乘估计（刚体估计职责）
  icp.h      点到点 ICP 主入口（迭代控制职责）
src/
  kdtree.cpp / rigid.cpp / types.cpp / icp.cpp
examples/example.cpp   已知变换 + 离群点的调用示例
tests/test_icp.cpp     单元/集成测试（28 项检查，无第三方测试框架）
Makefile
```

## 构建与运行

环境：GCC 11.4、GNU Make，Eigen 3.4.0 位于 `third_party/eigen3`。

```bash
make            # 构建静态库、示例 build/icp_example、测试 build/icp_tests
make test       # 运行全部测试
make example    # 运行配准示例
make clean
```

在其他项目中使用：编译时加 `-Iinclude -Ithird_party/eigen3`，链接 `-licp3d`。

## 公共接口

```cpp
#include "icp3d/icp.h"

icp3d::PointCloud source, target;                 // std::vector<Eigen::Vector3d>
Eigen::Matrix4d initial = Eigen::Matrix4d::Identity();

icp3d::IcpConfig config;
config.maxIterations             = 50;    // 最大迭代轮数（> 0）
config.maxCorrespondenceDistance = 1.0;   // 最大对应距离（与坐标同单位，> 0）
config.translationTolerance      = 1e-6;  // 平移增量收敛阈值（> 0）
config.rotationTolerance         = 1e-6;  // 旋转增量收敛阈值，弧度（> 0）
config.rankEpsilon               = 1e-10; // 旋转可解的奇异值判据（> 0）

icp3d::IcpResult result =
    icp3d::alignPointToPoint(source, target, initial, config);

// result.transform          最终（或失败时当前）源->目标 4x4 变换
// result.iterations         实际迭代轮数
// result.numCorrespondences 最终有效配对数
// result.rms                最终变换重新匹配的均方根距离
// result.converged / result.reason / result.message
```

## 算法说明

每一轮迭代：

1. 用当前变换 `T` 变换全部源点 `x' = R x + t`。
2. 在**只构建一次、各轮复用**的目标 KD 树上为每个变换后源点找最近目标点。
   允许多个源点匹配同一目标；等距时取**目标原始索引较小者**（结果确定）。
3. 剔除欧氏距离大于 `maxCorrespondenceDistance` 的配对；不保存完整距离矩阵，
   KD 树最近邻为 O(log n) 期望复杂度。
4. 对保留配对求最小二乘刚体增量 Δ（Horn/Umeyama，质心去均值 + SVD），
   反射情形翻转 `U` 的最后一列以禁止镜像，并累计 `T ← Δ·T`。
5. 平移增量 `‖Δt‖` 与旋转增量 `acos((tr(ΔR)−1)/2)` **同时**低于各自阈值才收敛。
6. 输出指标（配对数、RMS）一律用**最终变换重新匹配**计算，不沿用上一轮误差。

### 失败语义（不任取变换冒充成功）

| 场景 | 行为 |
| --- | --- |
| 点云为空 / 含非有限坐标 / 配置非正 / 初值非刚体 | 抛 `std::invalid_argument` |
| 有效配对少于 3 组 | `converged=false`，原因 `InsufficientCorrespondences` |
| 配对共线/共面等导致旋转不可唯一确定 | `converged=false`，原因 `DegenerateGeometry` |
| 达到最大迭代次数 | `converged=false`，原因 `MaxIterationsReached`，返回当前估计 |

### 数值容差

- 初始旋转校验：`max|RᵀR − I| ≤ 1e-9` 且 `|det R − 1| ≤ 1e-9`。
- 旋转可解性：去均值协方差矩阵的最小奇异值必须大于
  `rankEpsilon · max(1, σ₁)`（默认 `1e-10`）。平面/共线分布因此被判退化；
  这是合理的：绕分布平面法向等方向的旋转在数学上不可观。若点云本身为
  平面扫描，点到点 ICP 无法单独确定全部 3 个旋转自由度。
- 收敛判据使用严格小于阈值；角度经 `acos` 前夹紧到 `[-1, 1]`。
- 输出旋转由 SVD 正交矩阵乘积给出，天然满足正交与行列式 1，不引入缩放。

## 示例输出摘要

`make example` 使用 80 个已知变换采样点 + 2 个远距离离群点、240 个目标点：

```
配准前最近点 RMS: 9.686369
终止原因: Converged（已收敛）
迭代次数: 30, 最终有效配对数: 80（离群点已被距离阈值剔除）
配准后最近点 RMS（最终变换重新匹配）: 6.54e-15
```

## 测试覆盖

`tests/test_icp.cpp` 覆盖：KD 树正确性与等距 tie-break、SVD 刚体恢复、
正交/行列式保证、配对不足与共线退化、增量累计、镜像/缩放拒绝、ICP 收敛、
离群点剔除、配对不足失败、迭代上限未收敛、空集/非有限值/非法参数校验、
输入不可变。
