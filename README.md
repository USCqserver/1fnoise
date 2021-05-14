# 1/f noise

Simulation codes for quantum 1/f noise process in annealing and superconducting circuits. The theory is in Chapter 2 of my [thesis] (https://github.com/USCqserver/1fnoise/blob/master/Kawa_Yip_thesis.pdf). Codes included are:

:white_check_mark: spin_vector_model

:white_check_mark: qubit_model

:white_check_mark: CSFQ_circuit 

:white_check_mark: transmon_circuit (with dynamical decoupling)

:white_check_mark: ibm_crosstalk


### Input parameters:

- tf = total anneal time

- ntraj = number of trajectories

- nd = number of fluctuators per noise decade

- dec = number of noise decade

- bmean = mean of fluctuator strength

- bvariance = variance of fluctuator strength 

An annealing schedule of DWave NASA 2000Q is included. For CSFQ circuit, the default circuit parameter values (i_c, c_shunt, alpha, amp, ...) are in the scurve paper [[1]](#1).

## References
<a id="1">[1]</a> 
Khezri, M., Grover, J. A., Basham, J. I., Disseler, S. M., Chen, H., Novikov, S., ... & Lidar, D. A. (2021). Anneal-path correction in flux qubits. npj Quantum Information, 7(1), 1-8.

