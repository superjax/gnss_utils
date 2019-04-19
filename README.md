# GNSS Utils
[![Build Status](https://api.travis-ci.com/superjax/gnss_utils.svg?branch=master)](https://travis-ci.com/superjax/gnss_utils)

This is a bunch of utilities for handling Raw GNSS Data, as well as an example ROS node for calculating position, velocity and clock bias with iterated least-squares.

I made this so that I could feel confident that I was properly handling the data before I put this into a much larger project.  I would not recommend actually using the output of the ROS node for anything - there are _much_ better ways of doing GNSS positionin, but it's a good first-principles example.

## Building
The main library depends only on Eigen.  If you don't want the test suite or the ROS node, you should be able to build the library alone (or add it as a subdirectory of a larger project).  Building is the typical cmake - make workflow, or you can just put this in the `src` directory of a catkin workspace.

If catkin is detected, then the `point_positioning_demo` ROS node will be built.
If gtest is detected, then the test suite will be built.

I have abstracted as much as I can using C++ classes.  This results in (what I think is) really clean and understandable code.  Projects like RTKLIB are amazing, but the old-style C code is a little difficult for me to read.  My intention with releasing this package is that someone else out there will be able to get running faster by modifying my code, as opposed to starting with porting RTKLIB.

## THANKS:
I used a lot of code from [RTKLIB](https://github.com/tomojitakasu/RTKLIB) and [gnss-sdr](https://github.com/osqzss/gps-sdr-sim).  I mainly used these projects as references and copied a few of their functions into the `reference_algorithms.cpp` file which I test against.

There is a short LaTeX beamer presentation in the `doc` folder which contains a derivation of the algorithm I'm using for Iterated least-squares.  You should be able to build it with a standard LaTeX toolchain.

## TODO:
GLONASS Support






