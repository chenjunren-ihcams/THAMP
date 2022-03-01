# THAMP
This project is the implementation of THAMP algorithm presented in "A phenogenetic axis that modulates clinical manifestation and treatment outcome in myeloid neoplasms"

## 1. SYSTEM REQUIREMENTS

### All software dependencies and operating systems (including version numbers):
+ RStudio

+ MAC OS 12.0.1

### Software versions:
+ R version 4.0.5 (2021-03-31)

+ Python version 3.7.0 

### Required Python libraries:
+ keras == 2.5.0

+ tensorflow == 2.5.0

+ sklearn == 1.0.1

+ pandas == 1.2.4

+ numpy == 1.19.2

+ matplotlib == 2.2.3

+ joblib == 1.0.1

### Required R libraries:
+ readxl == 1.3.1

+ tidyverse == 1.3.1

+ ggplot2 == 3.3.5

+ spatstat == 2.2-0

+ sparr == 2.2-15

+ RColorBrewer == 1.1-2

+ ggrepel == 0.9.1

+ ggnewscale == 0.4.5

+ optparse = 1.7.1

+ showtext == 0.9-4

+ Cairo == 1.5-12.2

+ pROC == 1.18.0

### Any required non-standard hardware:
Not Applicable


## 2. INSTALLATION GUIDE
Not Applicable


## 3. INSTRUCTION

### Neural network architecture selection

#### *Run ANN_dim_reduce.py from "src/" folder:*

* Before running, delete the '#' at the beginning of line 25;

#### *Then run dimension_selection.R:*

* Before running, set the root path (path to "src/") at line 22.

* Quality metrics for two-dimensional embedding with varying neural network architectures are saved in folder  "../output/dim_sel".

### Main analysis

#### *Run ANN_dim_reduce.py from "src/" folder:*

* Before running, delete the entire "../output" folder and add a '#' at the beginning of line 25.

* Set neu_start=65 at line 37.

* Two-dimensional embedding of patients is saved in "../output/2d_full_data" folder (100 bootstrapping experiments). 

* Five-fold cross-validation result is saved as "../output/cv_stats/units_65_bs_0.csv". 

#### *Run THAMP-dirichlet-smooth-sort.R:*
* Before running, set the root path (path to "src/") at line 34.

* Lines 44-69 read in and preprocess data.

* Lines 74-171 project all myeloid patients onto the Pan-Myeloid Axis.

* Lines 175-526 perform Dirichlet clustering to generate Fig. 3A and Fig. 3B. 
(For  both mutation clustering and karyotype clustering) 

* Lines 529-850 define helper functions for calculating kernel densities of diseases, traits, and gene mutations.

* Lines 854-1111 generate Figs. 5B and 5C, saved in "../output/MYE/2d" folder.

* Lines 1114-1289 perform 100 bootstrapping experiments on the sorting of diseases, traits, and gene mutations; results are save in "../output/MYE/" folder.

* Lines 1295-1458 performed stability analysis; the results are stored in  "../output/MYE/mean_test" folder.

* Lines 1462-1678 generate Figs. 5D and 5E.

* Lines 1683-1888 calculate correlation coefficients. 

* Use "../output/MYE/cor/cor_kar_gene_total.csv" to generate Fig. 3C.

* Use "../output/MYE/cor/cor_gene_trait.csv" to generate Fig. 6A.

* Use "../output/MYE/cor/AML_cor_gene_trait.csv" to generate Fig. 6B.

* Use "../output/MYE/cor/AML_cor_gene_trait.csv" to generate Fig. 6C.

* The gene groups used for computing Fig. 7 are from the analysis explained in Fig. 6A. 

* Lines 1896-2005 analyzes treatment outcome in AML. 

* Lines 1952-1960 generate Fig. 7B.

* Lines 1963-1971 generate Fig. 7D.

* Lines 1976-2002 generate Fig. 7E.

* Lines 2006-2057 analyzes treatment outcome in MDS. 

* Lines 2034-2038 generate Fig. 8A.

* Lines 2040-2048 generate Fig. 8B.

* Lines 2052-2057 generate Fig. 8C.


## 4. POSTSCRIPT
* Total expected run time on a 'normal' desktop computer: 8 hours
