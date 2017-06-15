# SAClib

An efficient C++ implementation of sequential action control (SAC) for rapid online control of nonlinear dynamic systems.  This header-only library computes constrained optimal controls online to drive nonlinear robots according to an objective function, i.e., it is similar to nonlinear receding horizon control.  

The code closesly follows the reference in the **Citations** section.

## Overview

**Requirements:** Boost.

**Documentation:** A pdf reference manual is included in the ./doc folder.  Use doxygen to update the project's documentation.

**Getting Started:** See the examples in the ./examples folder.  For instance, the ./examples/Cart_and_Pend_full/ folder provides an example which uses the SAC library to invert a Cart and Pendulum system.  

The 'user' subfolder in each example, e.g., ./examples/Cart_and_Pend_full/user, is where the user should implement (or override) any class definitions required to implement dynamics, linearizations, or specialty cost functions.  In the Cart_and_Pend_full example, the dynamics and linearizations of the cart pendulum are provided and SAC optimizes a default quadratic cost functional to solve for controls.

## Citations

*This code was created while working on the sequential action control method with Prof. Todd Murphey in his robotics laboratory at Northwestern University.*

The SAC library is based on:

```
Ansari, Alexander R., and Todd D. Murphey. "Sequential action control: closed-form optimal control for nonlinear and nonsmooth systems." IEEE Transactions on Robotics 32.5 (2016): 1196-1214.
```

*Please include a reference* if you use this software.  For convenience, the bibtex is provided:

```
@Misc{SACSoftware2017,
author =   {Alexander R Ansari and Todd D Murphey},
title =    {{SAClib}: {A} {C++} header-only library implementation of sequential action control {(SAC)}},
howpublished = {\url{http://github.com/alexansari101/saclib}},
year = {2017}
}
```
```
@article{sacTRO2016,
  title={Sequential action control: closed-form optimal control for nonlinear and nonsmooth systems},
  author={Alexander R Ansari and Todd D Murphey},
  journal={IEEE Transactions on Robotics},
  volume={32},
  number={5},
  pages={1196--1214},
  year={2016},
  publisher={IEEE}
}
```
