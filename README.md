# SSDi-Calculator
 
---------------

## Overview

**SSDi-Calculator** is a python script that can be used to calculate 1) a sexual size dimorphism index (SSDi) based on average male and female sizes, 2) a SSDi based on the average of all male-female pairwise comparisons, and 3) conduct a permutation test to determine if the species is significantly size-dimorphic. The input files can be either comma-separated (csv) or tab-delimited text files (see below). 

## Version

The current release of **SSDi-Calculator** is v1.0.

## Installation

**SSDi-Calculator** is a module written in Python (3.7, and compatible with 2.7) that functions as a stand-alone command-line script (`SSDi-Calculator.py`). It can be downloaded and executed independently without the need to install **SSDi-Calculator** as a Python package or library, making it easy to use and edit. There are two Python packages that must be installed prior to use of **SSDi-Calculator**:

+ [numpy](http://www.numpy.org/)
+ [scipy](https://www.scipy.org/)

These packages (if not already shipped with your python) can be installed using common managers such as `pip` or `brew`. You can test that numpy and scipy are installed correctly by first opening python on the command line, and then importing them:

```
$ python
Python 3.7.3 (default, Jun 19 2019, 07:38:49) 
[Clang 10.0.1 (clang-1001.0.46.4)] on darwin
Type "help", "copyright", "credits" or "license" for more information.
>>> import numpy
>>> import scipy
>>> 
>>> quit()
```

If they are installed correctly, nothing obvious will happen on the screen. If they are not installed correctly, you will see an error:

```
>>> import Bio
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
ImportError: No module named numpy
```


## Citation

**SSDi-Calculator** is described the following publication:

+ Portik, D.M., Blackburn, D.C., and J.A. McGuire. 2020. Macroevolutionary patterns of sexual size dimorphism among African tree frogs (Family: Hyperoliidae). Journal of Heredity, Online Early. https://doi.org/10.1093/jhered/esaa019. 

[**Click here for PDF.**](https://academic.oup.com/jhered/advance-article/doi/10.1093/jhered/esaa019/5864412?guestAccessKey=25c9888e-e22c-4245-b78f-1722a37bcbbc)

If you use **SSDi-Calculator** for your research, please cite this publication.


## Instructions

### Quick Navigation:

+ [**1. Analysis Overview**](#AO)
+ [**2. Input Data Format**](#IDF)
+ [**3. Running the Script**](#RTS)
+ [**4. Example Data Set**](#EDS)

### 1. Analysis Overview <a name="AO"></a>

#### Calculating SSDi

This analysis will calculate a sexual size dimorphism index (SSDi) to quantify SSD across all species in the input file. The SSDi is calculated using the equation `SSDi = [(larger sex / smaller sex) – 1]`, arbitrarily set negative if males are the larger sex and positive if females are the larger sex (Lovich & Gibbons 1992). This SSDi has been widely used, is properly scaled around zero, and has high intuitive value because positive SSDi values indicate female-biased SSD and negative values indicate male-biased SSD. 

Two versions of the SSDi are calculated (for each species):

1. **Standard SSDi.** The mean body size of each sex is obtained and the means are used to calculate SSDi. We can call this the standard SSDi, as this is the typical way of calculating an SSDi. If you only have one data point for each sex for a given species, this is also fine. The point estimates will be used instead of the average.

2. **Pairwise SSDi.** The SSDi is obtained for all possible male-female pairwise comparisons, and the average SSDi is obtained from this distribution of values. We can call this the pairwise SSDi, as it results from pairwise comparisons. This version explicitly accounts for intraspecific variation. **If you only have one data point for each sex, this calculation is skipped.**  

#### Permutation Test

A permutation test is conducted if the species has data available for >1 male or >1 female. **If you only have one data point for each sex, this analysis is skipped.**  

When SSDi is calculated, a species is automatically classified as displaying male or female-biased SSD, even if the difference between male and female body size is minimal. In other words, a statistical distinction between sexual size dimorphism and monomorphism is not made based on this measure. The permutation test allows us to address this.

The permutation test is conducted with 10,000 bootstrap replicates to evaluate the null hypothesis that male and female body sizes come from the same population (where SSDi = 0). For each replicate:

1. The labels of males and females are randomly shuffled.

2. All possible pairwise SSDi values are calculated. 

3. A mean pairwise SSDi value is obtained from the set of values. 

The 10,000 simulated SSDi values of the permuted data represent the estimate of the sampling distribution under the null hypothesis. The test then assesses whether the empirical mean pairwise SSDi is outside of the critical values of the simulated distribution (2.5% and 97.5%, representing a 5% significance level), and a P-value is calculated. 

**If the P-value < 0.05:**
+ We can reject the null hypothesis that the sexes are equal in size, and classify the species as displaying sexual size dimorphism. 

**If the P-value > 0.05:**
+ We fail to reject the null hypothesis that the sexes are equal in size, and cannot classify the species as sexually size dimorphic. Although it would be tempting to say the species displays sexual size monomorphism, it is simply a convenience label. The failure to reject the null hypothesis could result from a true lack of body size difference between the sexes, or from artifacts such as insufficient sample sizes or high variation. Therefore, use caution when discussing these species and provide the caveats above.

### 2. Input Data Format <a name="IDF"></a>

The input can be either a CSV file (comma-separated) or tab-delimited text file. There are only three required columns:

1. **Species** - The species name. This can appear as a binomial name (`Afrixalus dorsalis` or `Afrixalus_dorsalis`, note spaces or undescores can be used) or a trinomial name (`Afrixalus dorsalis dorsalis`). Please note that names are not checked for errors (such as extra spaces), they are simply read as is. This means any differences in names (even small ones) will result in them being treated as different species. You can also use numbers in this column, as I did for the publication (in which lineages are denoted by integers - `Afrixalus dorsalis 1`, `Afrixalus dorsalis 2`, etc.).

2. **Sex** - This represents whether the individual is male or female. This should be coded as `m`, `M`, `f`, or `F`. The `m` and `f` are not case-sensitive (`m` = `M`). Any letters or symbols here that are not `m` or `f` will be excluded from the analysis (and written to the log file as errors).

3. **Size** - This should be the body size measurement of the individual, such as `35.5` or `35`. You can include any number of decimal places, or none at all.

Any number of additional columns can be included after the third column - they are simply ignored. The header line is also ignored, so you can label your columns however you like.

Here is an example of a tab-delimited input file:

```
Species	Sex	SUL(mm)	Museum Number	Source	Data Type	Country
Acanthixalus sonjae	f	35.5		Rodel et al. 2003	average	
Acanthixalus sonjae	m	35.1		Rodel et al. 2003	average	
Acanthixalus spinosus	m	33.1	CAS 153799	Present Study	specimen	Cameroon
Acanthixalus spinosus	f	33.7	CAS 153800	Present Study	specimen	Cameroon
Acanthixalus spinosus	f	34.5		Rodel et al. 2003	average	
Acanthixalus spinosus	m	33.5		Rodel et al. 2003	average	
Afrixalus brachycnemis	f	27.0		Channing 2001	maximum	
Afrixalus brachycnemis	m	25.0		Channing 2001	maximum	
Afrixalus dorsalis 1	m	24.5	CAS 253854	Present Study	specimen	Cameroon
Afrixalus dorsalis 1	f	26.9	CAS 253855	Present Study	specimen	Cameroon
Afrixalus dorsalis 1	m	21.0	CAS 253856	Present Study	specimen	Cameroon
Afrixalus dorsalis 1	m	24.4	CAS 253857	Present Study	specimen	Cameroon
Afrixalus dorsalis 1	m	24.0	CAS 253858	Present Study	specimen	Cameroon
Afrixalus dorsalis 1	m	26.3	CAS 253859	Present Study	specimen	Cameroon
Afrixalus dorsalis 1	m	26.0	CAS 253860	Present Study	specimen	Cameroon
```

And here is the same file but in CSV format:

```
Species,Sex,SUL (mm),Museum Number,Source,Data Type,Country
Acanthixalus sonjae,f,35.5,,Rodel et al. 2003,average,
Acanthixalus sonjae,m,35.1,,Rodel et al. 2003,average,
Acanthixalus spinosus,m,33.1,CAS 153799,Present Study,specimen,Cameroon
Acanthixalus spinosus,f,33.7,CAS 153800,Present Study,specimen,Cameroon
Acanthixalus spinosus,f,34.5,,Rodel et al. 2003,average,
Acanthixalus spinosus,m,33.5,,Rodel et al. 2003,average,
Afrixalus brachycnemis,f,27.0,,Channing 2001,maximum,
Afrixalus brachycnemis,m,25.0,,Channing 2001,maximum,
Afrixalus dorsalis 1,m,24.5,CAS 253854,Present Study,specimen,Cameroon
Afrixalus dorsalis 1,f,26.9,CAS 253855,Present Study,specimen,Cameroon
Afrixalus dorsalis 1,m,21.0,CAS 253856,Present Study,specimen,Cameroon
Afrixalus dorsalis 1,m,24.4,CAS 253857,Present Study,specimen,Cameroon
Afrixalus dorsalis 1,m,24.0,CAS 253858,Present Study,specimen,Cameroon
Afrixalus dorsalis 1,m,26.3,CAS 253859,Present Study,specimen,Cameroon
Afrixalus dorsalis 1,m,26.0,CAS 253860,Present Study,specimen,Cameroon
```

Notice that the column names can be whatever you'd like, and you can include as many columns after the first three. However, the order of the first three columns MUST be: species, sex, size.

Example files of both formats are provided in the [example-data](https://github.com/dportik/SSDi-Calculator/tree/master/example-data) folder.


### 3. Running the Script <a name="RTS"></a>

The script has three arguments that need to be supplied:

```
python SSDi-Calculator.py -i <input file> -f <input file format> -o <output directory>
```

#### Argument Explanations:

##### `-i <path-to-file>`

> The full path to a text file containing the data.

##### `-f <tab or csv>`

> The format of the input text file, either tab-delimited (tab) or comma-separated values (csv). You must enter one of these two choices: tab, csv.

##### `-o <path-to-directory>`

> The full path to an existing directory to write the output files.


#### Example Usage:

```
python SSDi-Calculator.py -i bin/Analysis/Hyperoliid-Dataset.csv -f csv -o bin/Analysis/Output/
```

This will treat `Hyperoliid-Dataset.csv` as a csv file, run the analyses, and write the output files to the `bin/Analysis/Output/` directory.

#### Output Files:

There are two output files produced.

**Results File**

The `SSDi-Results.csv` or `SSDi-Results.txt` files will contain all the results of the analysis. These files will contain the following column names:

+ `Species` - The name of the species.
+ `Number_Males` - The number of males with available size data.
+ `Number_Females` - The number of females with available size data.
+ `Avg_Male` - The mean body size of males.
+ `Avg_Female` - The mean body size of females.
+ `Standard_SSDi` - The SSDi calculated from the mean male and mean female size.
+ `Avg_Pairwise_SSDi` - The average SSDi value obtained from all pairwise SSDi values.
+ `AbsDifference` - The absolute difference between the standard SSDi and pairwise SSDi.
+ `Dimorphism_PValue` - The P-value associated with the permutation test.
+ `2.5_percentile` - The lower critical value obtained in the permutation test.
+ `97.5_percentile` - The upper critical value obtained in the permutation test.

Here is an example of the contents:

```
Species	Number_Males	Number_Females	Avg_Male	Avg_Female	Standard_SSDi	Avg_Pairwise_SSDi	AbsDifference	Dimorphism_PValue	2.5_percentile	97.5_percentile
Acanthixalus sonjae	1	1	NA	NA	0.011	NA	NA	NA	NA	NA
Acanthixalus spinosus	2	2	33.3	34.1	0.024	0.024	0.0	<0.001	-0.024	0.024
Afrixalus brachycnemis	1	1	NA	NA	0.08	NA	NA	NA	NA	NA
Afrixalus dorsalis 1	25	5	25.6	25.9	0.012	0.005	0.007	0.9	-0.084	0.061
Afrixalus dorsalis 2	10	2	24.4	26.4	0.082	0.087	0.005	0.123	-0.101	0.117
Afrixalus dorsalis 3	23	7	25.7	27.3	0.062	0.062	0.0	0.013	-0.052	0.046
Afrixalus dorsimaculatus 1	29	6	23.1	25.4	0.1	0.104	0.004	<0.001	-0.054	0.063
Afrixalus fornasini	32	11	31.0	32.7	0.055	0.054	0.001	<0.001	-0.025	0.026
Afrixalus fulvovittatus 1	35	8	24.4	25.9	0.061	0.068	0.007	0.076	-0.078	0.073
Afrixalus fulvovittatus 2	7	6	23.1	24.6	0.065	0.063	0.002	0.011	-0.051	0.051
```

Here you can see that several species contain only 1 data point for males and 1 data point for females - in these cases the pairwise SSDi and permutation tests are skipped and a `NA` is written to these columns. It also will display an `NA` for the average size columns, since only 1 data point is available for each and a mean cannot be calculated.

**Log File**

There will also be a log file produced called `SSDi-Calculator-Run.log`, with the same information printed to screen (and more!). This can be extremely useful for de-bugging and finding possible errors in your dataset. Here is an example of the contents:

```
DEBUG: 04-Aug-20 19:59:22: BEGINNING ANALYSIS
DEBUG: 04-Aug-20 19:59:22: Arguments: 
		-i /Users/portik/Documents/GitHub/SSDi-Calculator/example-data/Hyperoliid-Dataset.csv 
		-f csv 
		-o /Users/portik/Documents/GitHub/SSDi-Calculator/example-results
INFO: 04-Aug-20 19:59:22: Reading input data file: Hyperoliid-Dataset.csv
DEBUG: 04-Aug-20 19:59:22: input_to_dict: Skipping entry - Hyperolius concolor concolor, F?, 25.7
DEBUG: 04-Aug-20 19:59:22: input_to_dict: Skipping entry - Hyperolius concolor concolor, F?, 29.9
DEBUG: 04-Aug-20 19:59:22: input_to_dict: Skipping entry - Hyperolius fusciventris ssp, , 23.0
DEBUG: 04-Aug-20 19:59:22: input_to_dict: Skipping entry - Hyperolius fusciventris ssp, , 23.8
DEBUG: 04-Aug-20 19:59:22: input_to_dict: Skipping entry - Hyperolius guttulatus 2, F?, 33.9
DEBUG: 04-Aug-20 19:59:22: input_to_dict: Skipping entry - Hyperolius guttulatus 2, F?, 30.4
DEBUG: 04-Aug-20 19:59:22: input_to_dict: Skipping entry - Hyperolius guttulatus 2, F?, 27.6
DEBUG: 04-Aug-20 19:59:22: input_to_dict: Skipping entry - Hyperolius thomensis
DEBUG: 04-Aug-20 19:59:22: input_to_dict: Skipping entry - Hyperolius thomensis, F?, 40.4
DEBUG: 04-Aug-20 19:59:22: input_to_dict: Skipping entry - Tachycnemis seychellensis, J?, 39.9
INFO: 04-Aug-20 19:59:22: Found data for 164 species.


DEBUG: 04-Aug-20 19:59:22: Species with > 1 M and > 1 F: 95
DEBUG: 04-Aug-20 19:59:22: Species with > 2 M and > 2 F: 78
DEBUG: 04-Aug-20 19:59:22: Species with > 3 M and > 3 F: 50


INFO: 04-Aug-20 19:59:22: Species: Acanthixalus sonjae
INFO: 04-Aug-20 19:59:22: Males: 1
INFO: 04-Aug-20 19:59:22: Females: 1
INFO: 04-Aug-20 19:59:22: Standard SSDi: 0.011

INFO: 04-Aug-20 19:59:22: Species: Acanthixalus spinosus
INFO: 04-Aug-20 19:59:22: Males: 2
INFO: 04-Aug-20 19:59:22: Females: 2
INFO: 04-Aug-20 19:59:22: Standard SSDi: 0.024.
INFO: 04-Aug-20 19:59:22: Pairwise Analyses: Average pairwise SSDi: 0.024.
DEBUG: 04-Aug-20 19:59:22: Pairwise Analyses: One-sample t-test P-value: 0.053.
INFO: 04-Aug-20 19:59:22: Permutation Test: 2.5 and 97.5 percentile values: -0.024, 0.024.
INFO: 04-Aug-20 19:59:22: Permutation Test: Empirical value: 0.024
INFO: 04-Aug-20 19:59:22: Permutation Test: Empirical value lies outside the 97.5 percentile.
INFO: 04-Aug-20 19:59:22: Permutation Test: Permutation test P-value: <0.001.
```

Notice that some entries were skipped due to illegal values for the sex. The log file will also have results for a one-sample T-test (where the distribution of pairwise SSDi values is compared to a hypothesized mean of 0). The reviewers rightly pointed out this test will contain an inflated number of degrees of freedom, thereby skewing the results. That's why we use the permutation test instead. 

### 4. Example Data Set <a name="EDS"></a>

There is an example data set provided in the [example-data](https://github.com/dportik/SSDi-Calculator/tree/master/example-data) folder. The dataset is provided in tab-delimited format, and in CSV format. You can run **SSDi-Calculator** on these test files to ensure everything is working properly. The outputs from these datasets are provided in the [example-results](https://github.com/dportik/SSDi-Calculator/tree/master/example-results) folder.



## License

GNU Lesser General Public License v3.0

## Contact

SSDi-Calculator is written and maintained by Daniel Portik (daniel.portik@gmail.com). If you have any questions or issues, please open a topic on the `Issues` section on this github page.

