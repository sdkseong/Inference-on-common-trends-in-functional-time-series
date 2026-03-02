---
title: "Replication Code: Inference on Common Trends in Functional Time Series"
output: html_document
---

# Overview

[cite_start]This document outlines the replication procedures for the simulations and empirical analysis presented in the paper. [cite: 1]

---

## 1. Simulation Replication

[cite_start]The following scripts are used to replicate the statistical performance tables: [cite: 1]

### Size and Power (Table 1)
* [cite_start]**Execution:** Run `NSS_VRCode_Sim_SizePower.R` to replicate the size and power results. [cite: 1]
* [cite_start]**Reporting:** Run `report_table_SizePower.R` to generate the formatted results. [cite: 1]

### Dimension Estimation (Table 2)
* [cite_start]**Execution:** Run `NSS_VR_CODE_Sim_DIMEST.R` to replicate the correct dimension estimation results. [cite: 1]
* [cite_start]**Reporting:** Run `Report_table_DimEst1.R` to generate the formatted results. [cite: 1]

---

## 2. Yield Curve Empirical Analysis (Tables 3 & 4)

[cite_start]The analysis of yield curves utilizes raw data available at [https://home.treasury.gov](https://home.treasury.gov). [cite: 1]

### Model Generation
[cite_start]Run `NSS_VR_YieldCurve_main.R` to replicate Table 3 and Table 4. [cite: 1] [cite_start]The final commands in this script generate seven distinct models based on the following components: [cite: 1]

* [cite_start]$\varsigma_0$: **Level** [cite: 1, 2]
* [cite_start]$\varsigma_1$: **Slope** [cite: 1, 2]
* [cite_start]$\varsigma_2$: **Curvature** [cite: 1, 2]

### Model Descriptions

| Model | Null Hypothesis ($H_0$) | Component Representation |
| :--- | :--- | :--- |
| **Model 1** | Baseline | [cite_start]Reference model [cite: 1] |
| **Model 2** | $sp\{\varsigma_0\}$ | [cite_start]Level [cite: 1] |
| **Model 3** | $sp\{\varsigma_0, \varsigma_1\}$ | [cite_start]Level and Slope [cite: 1] |
| **Model 4** | $sp\{\varsigma_0, \varsigma_1, \varsigma_2\}$ | [cite_start]Level, Slope, and Curvature (Corresponds to **Table 4**) [cite: 1] |
| **Model 5** | $sp\{\varsigma_1\}$ | [cite_start]Slope [cite: 1] |
| **Model 6** | $sp\{\varsigma_2\}$ | [cite_start]Curvature [cite: 1] |
| **Model 7** | $sp\{\varsigma_1, \varsigma_2\}$ | [cite_start]Slope and Curvature [cite: 1] |

---
