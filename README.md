The current folder includes R code for reproducing all of the tables and figures in the article “Sample size requirements to test subgroup average treatment effects in cluster randomized trials” by Wang et al.

For questions or comments about the code please contact Xueqi Wang at xueqi.wang@yale.edu.

I. List of Supporting Files: These supporting files are sourced in the main files that reproduce the numbers in the submitted manuscript.
1. gendata.R = function to generate the simulated data sets for linear mixed modelling;
2. calcSubgroup.R = two functions to calculate the required number of clusters for the omnibus test and the intersection-union test respectively;
3. powerSubgroup.R = two functions to calculate the power for the omnibus test and the intersection-union test respectively.

II. List of Main Files: These main files are used to reproduce the results in the submitted manuscript.
4. empOmnibus.R = reproduce empirical type I error rate, empirical power, and predicted power results for the omnibus test in Table 1;
5. empIU.R = reproduce empirical type I error rate, empirical power, and predicted power results for the intersection-union test in Table 2;
6. conRES.R = reproduce sample size and power results for a back-of-the-envelope approach in Tables 1-2;
7. numStudy.R = reproduce results in Figures 1-2;
8. UMDEX.R = reproduce results of the application in Figures 3-4.

III. Folder
9. conResults = folder to save power/size results. Now this folder also contains 2 XLSX files of power/size results that could be reproduced through empOmnibus.R and empIU.R.

IV. Software
Analyses were conducted with R, version 4.2.2 (https://www.r-project.org/). The calculations used R packages nlme (version 3.1-160), openxlsx (version 4.2.5.1), mvtnorm (version 1.1-3), ggplot2 (version 3.4.0), gridExtra (version 2.3) directlabels (version 2021.1.13), and cowplot (version 1.1.1).

V. R commands for the installation of R packages
install.packages(c("nlme ", "openxlsx ", "mvtnorm ", "ggplot2", "gridExtra", "directlabels", "cowplot"))

NOTE: Make sure the current working directory is the current folder before running each program.
