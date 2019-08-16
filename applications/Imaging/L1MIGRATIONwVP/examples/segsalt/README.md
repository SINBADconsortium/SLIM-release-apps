The author recommends that all codes be run in parallel. Using 16 workers with 16GB of RAM per worker is recommended. If you run the codes interactively, remember to start matlabpool/parpool.

Approximate runtime on 8 local workers using Intel(R) Xeon(R) CPU E5-2680 v2 @ 2.80GHz:

* Fig2a_ideal_RTM_trueQ - under 1 hour
* Fig2b_ideal_Inv_trueQ - under 73 hours
* Fig2c_ideal_Inv_wrongQ - under 73 hours
* Fig2d_ideal_Inv_estQ - under 6 hours
* Fig4a_iwave_RTM_trueQ - under 1 hour
* Fig4b_iwave_Inv_trueQ - under 3 hours
* Fig4c_iwave_Inv_estQ - under 3 hours
