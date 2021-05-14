# 1/f noise

Simulation codes for quantum 1/f noise process, and its resulting dynamics with or without dynamical decoupling

:eight_spoked_asterisk: spin_vector_model

:eight_spoked_asterisk: qubit_model

:eight_spoked_asterisk: CSFQ_circuit 

:eight_spoked_asterisk: transmon_circuit 

:eight_spoked_asterisk: ibm_crosstalk


Input parameters:

:white_check_mark: tf = total anneal time

:white_check_mark: ntraj = number of trajectories

:white_check_mark: nd = number of fluctuators per noise decade

:white_check_mark: dec = number of noise decade

:white_check_mark: bmean = mean of fluctuator strength

:white_check_mark: bvariance = variance of fluctuator strength 

An annealing schedule of DWave NASA 2000Q is included. For CSFQ circuit, the default circuit parameter values (i_c, c_shunt, alpha, amp, ...) are in the scurve paper [[1]](#1).

## References
<a id="1">[1]</a> 
Khezri, M., Grover, J. A., Basham, J. I., Disseler, S. M., Chen, H., Novikov, S., ... & Lidar, D. A. (2021). Anneal-path correction in flux qubits. npj Quantum Information, 7(1), 1-8.

