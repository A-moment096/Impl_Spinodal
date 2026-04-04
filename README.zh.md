## 用不同编程语言实现调幅分解

探索相场模型在不同编程语言中实现可能性的爱好项目。

**其他语言版本：** [English](README.md) | [中文](README.zh.md)

### 介绍

本仓库包含A-B合金调幅分解（Spinodal Decomposition）的代码实现，使用了不同的编程语言编写。该仓库应与作者的个人博客 [A Moment's Rest](https://a-moment096.github.io) 同步更新。模拟案例来自 **S. Bulent Biner** 的《Programming Phase-Field Modeling》的第一个案例，但可能会在多个方面进行扩展。

本项目主要用于*存档博客文章中的代码*。关于代码的详细说明，请参考个人博客（目前仅提供中文版本）。

### 内容

目前包含以下内容：

- C++ 对应文章：[相场模拟，使用多种编程语言实现 I](https://a-moment096.github.io/posts/impl_spinodal_1/)
  - CPP_impl_v1.cpp: 最基础的实现，展示了一些 C++ 语言特性。
  - CPP_impl_v2.cpp: 更复杂的实现，包含面向对象特性和自定义边界条件。
- Python 对应文章：[相场模拟，使用多种编程语言实现 II](https://a-moment096.github.io/posts/impl_spinodal_2/)
  - PY_impl_v1.py: 使用基础 Python 数据结构的实现。
  - PY_impl_v2.py: 用 Numpy 数组替代 Python List，导致性能大幅下降。
  - PY_impl_v3.py: 使用 `numpy.roll` 函数实现周期性边界条件的高效算法。
  - PY_impl_v4.py: 添加 matplotlib 可视化计算结果，并使用 Python 面向对象特性。

### 未来计划

该系列计划涵盖以下编程语言（但不一定按顺序实现）。

- JavaScript
- TypeScript
- C
- Java
- Rust
- GoLang
- Haskell
- Bash Script
- ...

### 注意事项

以下是阅读（以及可能使用）本项目的некоторых注意事项。

1. 本项目纯粹是一个*爱好项目*，不保证代码的可用性。代码质量可能不高，因为它可能是在作者学习这些语言的过程中编写的。如果您发现任何bug，请谅解。
2. 本项目主要用于存档出现在个人博客中的代码。因此，除非有特殊原因（如严重的bug），现有代码不会被修改。
3. 不过，如果您有任何想法、问题或建议，欢迎分享！我们真诚欢迎 Issues 和 Pull Requests！
4. 如果您对本项目感兴趣并想使用它，请注意本项目采用 MIT License。感谢您的理解和支持。
