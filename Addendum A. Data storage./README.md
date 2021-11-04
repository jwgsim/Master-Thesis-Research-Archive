This Readme file contains:

1. Research description.
This folder contains the research manuscript which provides the reader with a 
very brief description of the hypotheses, research design, conceptual framework, method followed and 
structure of the data collected. 

2. Data files.
This folder contains two folders: A folder named "Empirical data" and a folder named 
"Simulation study data". 

The "Empirical data" folder contains the friendship relations and covariate information of the students in the 8th and 
16th classrooms of the Vermeij data. More specifically, this folder contains two files named "cov8.dat" and "cov16.dat",
which contain the covariate information for the students in the two classrooms. It additionally contains two files called 
"not8.net" and "not16.net", which contain the adjacency matrices describing the friendship relations between the students 
in the two classrooms. 

The "Simulation data" folder contains two R-data files, one called "simdata.RData" and one called "simdata_1_1.RData". 
The "simdata.RData" file contains data on 124 meta-analyses which have been obtained over each of the 12 simulation cells. 
More specifically, when loading the file into R, two RData files are loaded into the global environment. The first file 
contains information for the 124 meta-analyses over the three model complexities within the homogeneous and heterogeneous 
sample for the 25 sample size condition, where the second file contains information for the 124 meta-analyses over the 
three model complexities within the homogeneous and heterogeneous sample for the 75 sample size condition. The information 
that is provided for each seperate cell is:
- next counter: number of times the ERGM failed to converge. 
- model_GOF: goodness-of-fit indices for model statistics.
- Several _aux variables: goodness-of-fit indices for auxilary statistics.
- aic and bic: AIC and BIC for each ERGM fitted to one of the networks in the sample. 
- ma_model: meta-analysis results over the network sample.
- bias_matrix: matrix used to calculate bias and rmse.
- bias: bias for each meta-analysis parameter.
- rmse: rmse for each meta-analysis parameter.
- seed1: random seed used to estimate ERGM over each network in the sample. Can be used to reproduce presented results. 
- seed2: random seed used to obtain GOF over each network in the sample. Can be used to reproduce presented results.
The "simdata_1_1.RData" is finally required for generating plots 8 and 9 in the manuscript, because these contain the 
observed values in those two respective plots. More specifically, this data file contains the observed auxilary statistics
values for each of the networks in a single sample over which a meta-analysis has been obtained. To keep the scope of the
data manageable, this data has been saved for this one network sample, but not for all 124. 

3. Computer code.
This folder contains three folders: A folder named "Empirical computer code", a folder named "Simulation computer code", 
and a folder named "Results computer code". 

The "Empirical computer code" folder contains two files: an "Empirical data script.R" R-file and an "myEnvironment.RData"
file. The reader can use the "Empirical data script.R" file to obtain the "myEnvironment.RData" file. In short the 
"Empirical data script.R" uses the empirical data in the "Empirical data" folder to estimate the two true ERGMs which 
will be used to generate network samples. This information is stored in the "myEnvironment.RData" and will be required 
in the simulation computer code. 

The "Simulation computer code" folder contains one file, namely the "Simulation computer code.R" R-file. This file loads the 
"myEnvironment.RData" and uses it to simulate ERGMs based on the empirical networks, which it then re-estimates based on the
simulation conditions at hand. The script ultimately produces the "simdata.RData" file. The script used R version 4.0.3 
(2020-10-10), and was executed on a Windows 10 operating system.

The "Results computer code" folder finally contains one file, namely "Results computer code.R" which can be used to obtain
the results in the manuscript. More specifically, each result is given in the script in the order of the manuscript. Note
that all provided data files, namely the empirical data, the simulation data files "simdata.RData" and "simdata_1_1.RData",
and finally the "myEnvironment.RData" file from the "Empirical computer code.R" script, need to be loaded in the environment
to obtain all results. 

Additional information:
All documents and files in the research archive are saved by Jan-Willem Simons. The research thesis manuscript was 
submitted on 10/05/2021. The data was simulated over the 26-04-2021 to 02-05-2021 period. The empirical data on which
the simulation study is based was collected by dr. Lotte Vermeij between December 1999 and May 2000. 
 
