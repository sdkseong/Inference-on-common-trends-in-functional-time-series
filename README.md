---
title: "Replication Code: Inference on Common Trends in Functional Time Series"
output: html_document
---

# Overview

This document outlines the replication procedures for the simulations and empirical analysis presented in the paper. 

---

## 1. Simulation Replication

The following scripts are used to replicate the statistical performance tables:

### Size and Power (Table 1)
**Execution:** Run `NSS_VRCode_Sim_SizePower.R` to replicate the size and power results. 
**Reporting:** Run `report_table_SizePower.R` to generate the formatted results.

### Dimension Estimation (Table 2)
**Execution:** Run `NSS_VR_CODE_Sim_DIMEST.R` to replicate the correct dimension estimation results.
**Reporting:** Run `Report_table_DimEst1.R` to generate the formatted results. 

---

## 2. Yield Curve Empirical Analysis (Tables 3 & 4)

The analysis of yield curves utilizes raw data available at [https://home.treasury.gov](https://home.treasury.gov).

### Model Generation
Run `NSS_VR_YieldCurve_main.R` to replicate Table 3 and Table 4. The final commands in this script generate seven distinct models based on the following components

*  $\varsigma_0$: **Level** 
*  $\varsigma_1$: **Slope**  
*  $\varsigma_2$: **Curvature**  

### Model Descriptions

| Model | Null Hypothesis ($H_0$) | Component Representation |
| :--- | :--- | :--- |
| **Model 1** | Baseline |  Reference model   |
| **Model 2** | $sp\{\varsigma_0\}$ |  Level   |
| **Model 3** | $sp\{\varsigma_0, \varsigma_1\}$ |  Level and Slope   |
| **Model 4** | $sp\{\varsigma_0, \varsigma_1, \varsigma_2\}$ | Level, Slope, and Curvature (Corresponds to **Table 4**)   |
| **Model 5** | $sp\{\varsigma_1\}$ | Slope   |
| **Model 6** | $sp\{\varsigma_2\}$ |  Curvature   |
| **Model 7** | $sp\{\varsigma_1, \varsigma_2\}$ |  Slope and Curvature  |

---
