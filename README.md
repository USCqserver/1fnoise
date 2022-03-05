<h1 align="center" style="display: block; font-size: 2.5em; font-weight: bold; margin-block-start: 1em; margin-block-end: 1em;">
<br><br><strong> 1/f noise</strong>
</h1>

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Latest release](https://img.shields.io/github/v/release/USCqserver/1fnoise)](https://github.com/USCqserver/1fnoise/releases/tag/v2.0.0)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2022a-BLUE.svg)](https://www.mathworks.com/login?uri=%2Fdownloads%2Fprerelease)

---
## Table of contents[![](./docs/img/pin.svg)](#table-of-contents)
1. [Introduction](#introduction)
2. [Code description](#codedescription)
3. [Input parameters](#inputparameters)
4. [Reference](#reference)
---
## Introduction[![](./docs/img/pin.svg)](#introduction)
Simulation codes for time-dependent quantum 1/f noise process in annealing and superconducting circuits. The theory is in Chapter 2 of my thesis [[1]](#1). 

## Code description <a name="codedescription"></a>
Codes included are:

:white_check_mark: [spin_vector_model](https://github.com/USCqserver/1fnoise/tree/master/spin_vector_model)

:white_check_mark: [qubit_model](https://github.com/USCqserver/1fnoise/tree/master/qubit_model)

:white_check_mark: [CSFQ_circuit](https://github.com/USCqserver/1fnoise/tree/master/CSFQ_circuit) 

:white_check_mark: [transmon_circuit](https://github.com/USCqserver/1fnoise/tree/master/transmon_circuit) (with dynamical decoupling)

:white_check_mark: [ibm_crosstalk](https://github.com/USCqserver/1fnoise/tree/master/ibm_crosstalk)


## Input parameters <a name="inputparameters"></a>

- `tf`: total anneal time

- `ntraj`: number of trajectories

- `nd`: number of fluctuators per noise decade

- `dec`: number of noise decade

- `bmean`: mean of fluctuator strength

- `bvariance`: variance of fluctuator strength 

A schedule of DW2X annealing schedule is included. For CSFQ circuit, the default circuit parameter values (`i_c`, `c_shunt`, `alpha`, `amp`, ...) are in the scurve paper [[2]](#2). An AME quantum trajectories code of the CSFQ circuit is also included.

## Reference <a name="reference"></a>
<a id="1">[1]</a> 
Yip, K.W., 2021. [Open-system modeling of quantum annealing: theory and applications](https://arxiv.org/pdf/2107.07231.pdf). arXiv preprint arXiv:2107.07231.

<a id="2">[2]</a> 
Khezri, M., Grover, J. A., Basham, J. I., Disseler, S. M., Chen, H., Novikov, S., ... & Lidar, D. A. 2021. 
[Anneal-path correction in flux qubits](https://arxiv.org/pdf/2002.11217.pdf). npj Quantum Information, 7(1), 1-8.

