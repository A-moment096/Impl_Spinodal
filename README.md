## Implment Spinodal with Different Programming Language

A hobby project to explore the possiblity of phase field with programming languages.

**Available in:** [English](README.md) | [中文](README.zh.md)

### Introduction

This repository contains code for spinodal decomposition of A-B alloy, written in different programming languages. This repository should be updated with my personal blog [A Moment's Rest](https://a-moment096.github.io). The simulation case is from *Programming Phase-Field Modeling*, case study 1 by **S. Bulent Biner**, but may be extended in several means.

This project is mainly for *archiving the code in the blog posts*. For explanation of the code, please refer to my blog (but currently it's in Chinese only, sorry for inconvenience).

### Contents

Currently it contains:

- C++ corresponding post: [Phase Field Simulation, but with may languages I](https://a-moment096.github.io/posts/impl_spinodal_1/)
  - CPP_impl_v1.cpp: most basic implementation with some C++ language features demonstrated.
  - CPP_impl_v2.cpp: more complex implementation with OOP feature and customise boundary conditions.
- Python corresponding post: [Phase Field Simulation, but with may languages II](https://a-moment096.github.io/posts/impl_spinodal_2/)
  - PY_impl_v1.py: implementation with primitive Python data structures.
  - PY_impl_v2.py: use Numpy array to replace `List` of Python, causing dramatic performance decrease.
  - PY_impl_v3.py: use `numpy.roll` function to implement efficient algorithm for periodic boundary condition calculation.
  - PY_impl_v4.py: add matplotlib to visualize the result on the fly after calculation finished and use OOP features of Python.

### Furture Plan

For this series, we are going to cover the following programming languages (but may not implement in this order though).

- JavaScript
- TypeScript
- C
- Java
- Rust
- GoLang
- Haskell
- Base Script
- ...

### Notice

Here are some notice for reading (and possibly using) this project.

1. This project is purly a *hobby project* and the code availablity is not granted. The code could be in bad quality as it could be written while author learning these languages. So, if you find any bug, please bear with me.
2. This project is mainly for archive the code appears in my blog. So without any special reason (like serious bugs), the code existed would not be modified. 
3. Although, if you have any idea, question or suggestion, welcome to share with us! Issues and PR are warmly welcome!
4. If you are interested in this project and would like to use it, please notice that this project is under MIT License. Thanks for your kindly consideration.
 