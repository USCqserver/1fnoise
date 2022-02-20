<h1 align="center" style="display: block; font-size: 2.5em; font-weight: bold; margin-block-start: 1em; margin-block-end: 1em;">
<br><br><strong> 1/f noise</strong>
</h1>

![Latest release](https://img.shields.io/github/v/release/aregtech/areg-sdk?label=%20%F0%9F%93%A3%20Latest%20release&style=flat&logoColor=b0c0c0&labelColor=363D44)
<img src="https://img.shields.io/badge/MATLAB-R2022a-BLUE.svg" alt="MATLAB solution"/>

---
## Table of contents[![](./docs/img/pin.svg)](#table-of-contents)
1. [Introduction](#introduction)
2. [Code description](#codedescription)
3. [Input parameters](#inputparameters)
4. [Reference](#reference)
---
## Introduction[![](./docs/img/pin.svg)](#introduction)
Simulation codes for time-dependent quantum 1/f noise process in annealing and superconducting circuits.<br /> 

The theory is in Chapter 2 of my [thesis](https://github.com/USCqserver/1fnoise/blob/master/Kawa_Yip_thesis.pdf). Codes included are:

:white_check_mark: spin_vector_model

:white_check_mark: qubit_model

:white_check_mark: CSFQ_circuit 

:white_check_mark: transmon_circuit (with dynamical decoupling)

:white_check_mark: ibm_crosstalk


### Input parameters: <a name="inputparameters"></a>

- `tf`: total anneal time

- `ntraj`: number of trajectories

- `nd`: number of fluctuators per noise decade

- `dec`: number of noise decade

- `bmean`: mean of fluctuator strength

- `bvariance`: variance of fluctuator strength 

A schedule of DW2X annealing schedule is included. For CSFQ circuit, the default circuit parameter values (`i_c`, `c_shunt`, `alpha`, `amp`, ...) are in the scurve paper [[1]](#1). An AME quantum trajectories code of the CSFQ circuit is also included.

## Reference <a name="reference"></a>
<a id="1">[1]</a> 
Khezri, M., Grover, J. A., Basham, J. I., Disseler, S. M., Chen, H., Novikov, S., ... & Lidar, D. A. (2021). Anneal-path correction in flux qubits. npj Quantum Information, 7(1), 1-8.

