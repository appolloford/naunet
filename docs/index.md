# Welcome to Naunet

## Introduction

Naunet is a general-purpose astrochemistry tool modeling the chemical evolution in the universe. It can work with the network from KIDA or UMIST databases. KROME users can also easily bring their original network files here. 

The idea of Naunet is acting as a code-generator to generate C++ codes according to the input chemical network. These codes are prepared for the `CVODE` differential equation solver in `SUNDIALS` package and `Odeint` in boost library. For the multiple methods implemented in `CVODE` and `Odeint`, Naunet is capable of generating different codes for the methods. These codes can be compiled into libraries to be flexibly used.

## System requirements

Naunet requires Python 3.8+, CMake 3.18+, and a compiler supporting C++11. Depending on the solver, users also have to install SUNDIALS > 5.6.1 or 1.65.0 < Boost < 1.75.0. We recommend users to start from SUNDIALS 5.7.0, which has been tested in the Github environment. Only CVODE is used in the projects, so it should be possible to only install CVODE. If you are using KLU solver, SuiteSparse is required by SUNDIALS. Similarly, if you are using the cuSolverSp_batchQR solver, then CUDA is required by SUNDIALS.

## Installation

Naunet can be installed by pip:
```console
~ $ pip install naunet
```
> CAVEAT: Naunet is still under developing. This may not install the latest version.

To install the latest version, clone the directory from Github and install from the directory:
```console
~ $ git clone https://github.com/appolloford/naunet
~ $ cd naunet
~ $ pip install .
```

## Getting started

To use Naunet, you will have to prepare a network. Here we will start from a prepared example network. First, let's create a folder and create the example in it.
```console
~ $ mkdir naunet_test && cd naunet_test
~ $ naunet example --select=0
```
> If you call `naunet example`, it will show a list of all available examples. You can pick one of them as you like.
