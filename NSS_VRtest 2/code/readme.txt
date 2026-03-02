################################################################
############ Replication Code ##############
#### "Inference on common trends in functional time series" #### 
################################################################


Simulation replication is available:

NSS_VRCode_Sim_SizePower.R : Replicate Size and Power Table (Table 1)
report_table_SizePower.R: Report results


NSS_VR_CODE_Sim_DIMEST.R :Replicate Correct dimension estimation Table (Table 2)
Report_table_DimEst1.R: Report results



NSS_VR_YieldCurve_main.R : Replicate Table 3 and Table 4
The raw data is available at https://home.treasury.gov  
Note that the last commands generate seven models 

## Model 1: Baseline
## Model 2: H_0 = sp\{\varsigma_0}
## Model 3: H_0 = sp\{\varsigma_0, \varsigma_1}
## Model 4: H_0 = sp\{\varsigma_0, \varsigma_1, \varsigma_2}
## Model 5: H_0 = sp\{\varsigma_1}
## Model 6: H_0 = sp\{\varsigma_2}
## Model 7: H_0 = sp\{\varsigma_1, \varsigma_2}

\varsigma_0, \varsigma_1, and \varsigma2 represent the level, slope, and curvature components respectively. 

The testing results for Model 4 corresponds to Table 4.
