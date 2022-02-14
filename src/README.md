# README for code on THAMP data

## 1. SYSTEM REQUIREMENTS

### All software dependencies and operating systems (including version numbers):
RStudio
MAC OS 12.0.1

### Software versions:
R version 4.0.5 (2021-03-31)
Python version 3.7.0 

### Required Python libraries:
keras == 2.5.0
tensorflow == 2.5.0
sklearn == 1.0.1
pandas == 1.2.4
numpy == 1.19.2
matplotlib == 2.2.3
joblib == 1.0.1

### Required R libraries:
readxl == 1.3.1
tidyverse == 1.3.1
ggplot2 == 3.3.5
spatstat == 2.2-0
sparr == 2.2-15
RColorBrewer == 1.1-2
ggrepel == 0.9.1
ggnewscale == 0.4.5
optparse = 1.7.1
showtext == 0.9-4
Cairo == 1.5-12.2

### Any required non-standard hardware:
Not Applicable


## 2. INSTALLATION GUIDE
Not Applicable


## 3. INSTRUCTION

### Neural network architecture selection

Run ANN_dim_reduce.py from "src/" folder:
Before running, delete the '#' at the beginning of line 25;

Then run dimension_selection.R:
Before running, set the root path at line 22.

Quality metrics for two-dimensional embedding with varying neural network architectures are saved in folder  "../output/dim_sel".

### Main analysis

Run ANN_dim_reduce.py from "src/" folder:
Before running, delete the entire "../output" folder and add a '#' at the beginning of line 25.
Set neu_start=65 at line 37.
Two-dimensional embedding of patients is saved in "../output/2d_full_data" folder (100 bootstrapping experiments). 
Five-fold cross-validation result is saved as "../output/cv_stats/units_65_bs_0.csv". 

Run THAMP-dirichlet-smooth-sort.R:
Before running, set the root path at line 34.
Lines 44-69 read in and preprocess data.
Lines 74-171 project all myeloid patients onto the Pan-Myeloid Axis.

Lines 175-521 perform Dirichlet clustering to generate Fig. 2a and Fig. 2b. 
(For mutation clustering, set "mutation_or_karyotype <- mutation" at line 181.
 For karyotype clustering, "mutation_or_karyotype <- karyotype" at line 181.) 

Lines 524-845 define helper functions for calculating kernel densities of diseases, traits, and gene mutations.
Lines 849-1106 generate Figs. 3a and 3b, saved in "../output/MYE/2d" folder.

Lines 1109-1284 perform 100 bootstrapping experiments on the sorting of diseases, traits, and gene mutations; results are save in "../output/MYE/" folder.
Lines 1290-1453 performed stability analysis; the results are stored in  "../output/MYE/mean_test" folder.
Lines 1457-1673 generate Figs. 4a and 4c.

Lines 1678-1883 test Simpson's Paradox. 
Use "../output/MYE/group/AML group/AML_pval_log.csv" to generate Fig. 5a.
Use "../output/MYE/group/MDS group/MDS_pval_log.csv" to generate Fig. 5b.

The gene groups used for computing Fig. 6 are from the analysis explained in Fig. 4b. 
Lines 1891-1930 analyzes treatment outcome in AML. 
Lines 1924-1926 generate Fig. 6a.
Lines 1928-1930 generate Supplementary Fig. 3.  

Lines 1933-1969 analyzes treatment outcome in MDS. 
Lines 1962-1964 generate Fig. 6d.
Lines 1967-1969 generate Supplementary Fig. 4.

## 4. POSTSCRIPT
Total expected run time on a 'normal' desktop computer:
8 hours




