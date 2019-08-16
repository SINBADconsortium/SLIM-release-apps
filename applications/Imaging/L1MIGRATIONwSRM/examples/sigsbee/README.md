The author recommends that all codes be run in parallel. Using 31 workers with 16GB of RAM per worker is recommended. If you run the codes interactively, remember to start matlabpool/parpool.

Approximate runtime on 8 local workers using Intel(R) Xeon(R) CPU E5-2680 v2 @ 2.80GHz:

* Fig4c_RTM_Total - under 6 hours
* Fig4d_Inv_Total - under 24 hours
* Fig4e_Inv_Total_IgnoreMul - under 23 hours
* Fig5a_Inv_Mul - under 24 hour
* Fig5b_Inv_Total_EstQ - under 52 hours
* Fig6c_Inv_Total_AccuBgModel - under 1 minute
* Fig6d_Inv_Total_SlowBgModel - under 26 hours
* Fig6e_Inv_Total_AccuBgModel_NoSparsity - under 20 hours
* Fig6f_Inv_Total_SlowBgModel_NoSparsity - under 20 hourS
