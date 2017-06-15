# SAClib

An efficient C++ implementation of sequential action control (SAC) for rapid online control of nonlinear dynamic systems.  This header-only library computes constrained optimal controls online to drive nonlinear robots according to an objective function, i.e., it is similar to nonlinear receding horizon control.  

The code closesly follows the paper referenced in the **Citations** section.

## Overview

**Requirements:** Boost (the examples use several of the optional compiled Boost libraries).  This software has been compiled and tested on GCC 4 (Ubuntu 14.04) and GCC 5 (Ubuntu 16.04).

**Documentation:** A pdf reference manual is included in the ./doc folder.  Use doxygen to update the project's documentation.

**Getting Started:** Install a local copy of saclib to your computer by cloning this repository.  

For examples of how to use this library see the examples in the ./examples folder.  For instance, the ./examples/Cart_and_Pend_full/ folder provides an example that uses the SAC library to invert an acceleration controlled cart and pendulum system.  This example uses a default, quadratic (LQR-like) state tracking objective that measures the error between the current state (configuration + velocities) of the cart and penulum system and the origin (inverted equilibrium).  To minimize this objective, SAC computes constrained acceleration controls that move the cart laterally in order to invert the pendulum.  

The 'user' subfolder in each example, e.g., ./examples/Cart_and_Pend_full/user, is where the user should implement (or override) any class definitions required to implement dynamics, linearizations, or specialty cost functions.  In the Cart_and_Pend_full example, the dynamics and linearizations of the cart pendulum are provided and SAC optimizes a default quadratic cost functional to solve for controls as described in the Transactions on Robotics Paper (see **Citations**).

Note: The user can override any header file in the library by including a user modified copy in an example's 'user' folder.

Remember to __edit the Makefile__ in each example folder so that they link to the installed location of the SAClib library.  Build each example by calling "make" from within the folder pertaining to the appropriate example.

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
howpublished = {\url{https://github.com/alexansari101/saclib}},
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

## Contributing

We encourage contributions to this library.  The library source files can be found in the ./lib/src folder.  Use Doxygen style comments.  Test code with appropriate unit tests and ensure the examples compile/run.

Contributors should add their names to the AUTHORS file.
