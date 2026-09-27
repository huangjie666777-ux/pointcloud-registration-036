# 工件扫描点云配准库

初始项目仅包含Eigen3.4.0头文件与构建环境说明，尚无业务实现。

- C++17，GCC11.4.0，可直接使用g++。
- GNU Make4.3，可直接使用make。
- Eigen3.4.0位于third_party/eigen3，编译时使用-std=c++17和-Ithird_party/eigen3。
- Eigen来自Ubuntu22.04的libeigen3-dev=3.4.0-2ubuntu2包，版权与许可见third_party/eigen3/COPYRIGHT。
