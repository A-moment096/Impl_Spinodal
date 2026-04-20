## 用不同编程语言实现调幅分解

探索相场模型在不同编程语言中实现可能性的爱好项目。

**其他语言版本：** [English](README.md) | [中文](README.zh.md)

### 介绍

本仓库包含A-B合金调幅分解（Spinodal Decomposition）的代码实现，使用了不同的编程语言编写。该仓库应与作者的个人博客 [A Moment's Rest](https://a-moment096.github.io) 同步更新。模拟案例来自 **S. Bulent Biner** 的《Programming Phase-Field Modeling》的第一个案例，但可能会在多个方面进行扩展。

本项目主要用于*存档博客文章中的代码*。关于代码的详细说明，请参考个人博客（但是目前仅提供中文版本）。

### 内容

目前包含以下内容：

- C++ 
  > 对应博客：[相场模拟，但是用很多语言 I](https://a-moment096.github.io/posts/impl_spinodal_1/)
  - `.gitignore`: 让 Git 忽略 `build` 文件夹。
  - `CMakeLists.txt`: CMake 项目文件，用来构建下面的两个源文件。
  - `CPP_impl_v1.cpp`: 最基础的实现，展示了一些 C++ 语言特性。
  - `CPP_impl_v2.cpp`: 更复杂的实现，包含面向对象特性和自定义边界条件。

- Python 
  > 对应博客：[相场模拟，但是用很多语言 II](https://a-moment096.github.io/posts/impl_spinodal_2/)
  - `PY_impl_v1.py`: 使用基础 Python 数据结构的实现。
  - `PY_impl_v2.py`: 用 Numpy 数组替代 Python List，导致性能大幅下降。
  - `PY_impl_v3.py`: 使用 `numpy.roll` 函数实现周期性边界条件的高效算法。
  - `PY_impl_v4.py`: 添加 matplotlib 可视化计算结果，并使用 Python 面向对象特性。

- JavaScript & TypeScript 
  > 对应博客: [相场模拟，但是用很多语言 III](https://a-moment096.github.io/posts/impl_spinodal_3/)
  - `.gitignore`: 让 Git 忽略 `.vscode` 和 `node_modules` 两个文件夹。
  - `JS_impl_v1.js`: 使用原生 JavaScript 实现 `1-C++/CPP_impl_v1.cpp` 中的计算。 
  - `JS_impl_v2.js`: 在第一版的基础上加入面向对象和函数式写法。
  - `JS_impl_v3_config.js`: 允许从 `simu_config.json` 文件中读取并解析模拟参数。 
  - `JS_impl_v3_data_structures.js`: 在 `JS_impl_v2.js` 中定义的大部分数据结构。
  - `JS_impl_v3.js`: 第三版的主入口。
  - `packages.json`: JS 和 TS 的项目文件。
  - `pnpm-lock.yaml`: `pnpm` 生成的 Lockfile，用来同步项目依赖。
  - `simu_config.json`: 可以被 `JS_impl_v3.js` 和 TypeScript 版本读取的模拟设置文件。
  - `TS_impl_v1.ts`: 用 TypeScript 实现的 `JS_impl_v3.js`。
  - `TS_impl_v2.ts`: 使用 React 以及 `TS_impl_v2.html` 来在浏览器中渲染模拟结果。
  - `TS_impl_v2.html`: 用于渲染 `TS_impl_v2.ts` 的计算结果并显示在浏览器中。
  - `tsconfig.json`: TypeScript 的配置文件。

- C 
  > 当前对应的博客还没有写好……
  - `.vscode`: VS Code 相关的配置文件。用来在 Linux 和 Windows 上构建、调试和运行项目。
  - `include`: 外部库用到的头文件（目前只有 `fftw3`）.
  - `lib`: 这个项目用到的库。包含自实现的 `myfft` 和下载的 `fftw3`.
    > 你可以使用下文的 `C_my_fft.c` 来自己编译出这里的 `myfft` 库。在 Windows 下编译 `C_impl_fft_v3.c` 并运行时，请将 `libfftw3-3.dll` 放在和编译产物相同的文件夹内，确保能正常运行。
  - `scripts`: 用来自行使用或让 VS Code 进行编译/调试的编译脚本。
  - `C_impl_fft_v1.c`: 使用傅立叶谱法解 Cahn-Hilliard 方程实现模拟。使用了自实现的 FFT 算法以及显示解法。默认使用递归 FFT 算法并输出 VTK 文件。
  - `C_impl_fft_v1.c`: 使用自实现的 FFT 算法以及半隐式解法。默认使用迭代 FFT 算法，不输出 VTK 文件。
  - `C_impl_fft_v1.c`: 使用 `fftw3` 代替自实现算法。采用半隐式解法，默认不输出 VTK 文件。
  - `C_my_fft_complex.h`: 实现类似 `fftw_complex` 的复数结构.
  - `C_my_fft_test_main.c`: `C_my_fft.c` 和 `C_my_fft_complex.h` 的单元测试。
  - `C_my_fft.c`: 实现递归和迭代 FFT 算法。采用 Cooley–Tukey FFT 算法。
  - `C_my_fft.h`: 使用 `C_my_fft.c` 或预编译的 `libmyfft.lib` 时要用的头文件。
  - `create_directory.h`: 跨平台创建文件夹的函数实现。
  - `platform.h`: 平台相关的设置。 包含了测速函数以及 `create_directory.h`.


### 未来计划

该系列计划涵盖以下编程语言（但不一定按顺序实现）。

- JavaScript（*已实现！*）
- TypeScript（*已实现！*）
- C（*已实现，还没写好*）
- Java
- Rust
- GoLang
- Haskell
- Bash Script
- ...

### 注意事项

以下是阅读（以及可能使用）本项目的注意事项。

1. 本项目纯粹是一个*爱好项目*，不保证代码的可用性。代码质量可能不高，因为它可能是在作者学习这些语言的过程中编写的。如果您发现任何bug，请谅解。
2. 本项目主要用于存档出现在个人博客中的代码。因此，除非有特殊原因（如严重的bug），现有代码不会被修改。
3. 不过，如果您有任何想法、问题或建议，欢迎分享！我们真诚欢迎 Issues 和 Pull Requests！
4. 如果您对本项目感兴趣并想使用它，请注意本项目采用 MIT License。感谢您的理解和支持。
