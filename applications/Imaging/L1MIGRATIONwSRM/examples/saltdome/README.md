The author recommends that all codes be run in parallel. Using 15 workers with 16GB of RAM per worker is recommended. If you run the codes interactively, remember to start matlabpool/parpool.

Approximate runtime on 8 local workers using Intel(R) Xeon(R) CPU E5-2680 v2 @ 2.80GHz:

* Fig1c_RTM_IgnoreMul - under 1 hour
* Fig1d_RTM_totaldata - under 1 hour
* Fig2a_Inv_alldata - under 130 hours
* Fig2bc_Inv_RTM_10ss - under 18 hours
* Fig2d_Inv_2ss15freq_NoRedraw - under 12 hours
* Fig2e_Inv_2ss15freq - under 8 hours
* Fig2f_Inv_2ss15freq_L2P - under 8 hours
