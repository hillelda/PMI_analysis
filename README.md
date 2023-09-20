# PMI Analysis

Recent research suggests that Post-Mortem Interval (PMI) plays a critical role when sequencing mRNA.
This project aims to investigate whether the same holds true for tRNA fragments (tRFs), which have been implicated in 
regulating protein synthesis within cells. Although evidence for the relevance of PMI to tRF counts was found in liver tissue,
this project analyzed brain tissue, where much less enzymes are found, and therefore we assumed that it is different from liver tissue. 

### Objective

The primary objective of this analysis is to determine whether there is a correlation between PMI and the count of tRFs
in patient samples. Specifically, we are interested in understanding if higher PMI values result in greater degradation
of tRFs. If this hypothesis is correct, we would expect the density distribution of tRFs for different tRF groups
(3', 5', i, 5-half, 3-half) to shift to the left, indicating a higher density of shorter tRFs.

### Analysis Overview

Step-by-step breakdown of the analysis conducted:

    1. Data Gathering: The code gathers patient information, including their tRF sequencing data, to be used for analysis.

    2. Data Grouping: The tRF counts are grouped based on PMI and Alzheimer's Disease (AD) condition. This grouping helps in exploring the relationships between tRFs, PMI, and AD.

    3. Density Plot: The analysis includes the creation of density plots to visualize the distribution of tRFs within different patient groups. These plots provide insights into any potential correlations.

### Results

Upon examining the density plots, we did not observe any significant shifts that would indicate a correlation between tRFs and PMI. This suggests that, unlike mRNA, tRFs may not be significantly affected by PMI in the studied samples.

# Information 
### Dataset Information

The dataset used in this analysis was obtained from RUSH institution. It contains tRF sequencing data and metadata for 48 patients.
The dataset is available on request from the RUSH institution, and is not provided in this repository.

### Project Team

This analysis was conducted at Soreq's lab, Edmond and Lily Safra Center for Brain Sciences; under the guidance of Prof. Hermona Soreq.
The project team involved in this analysis includes Hillel Darshan, Nimrod Madrer, and  Dr. Daived Greenberg. 
