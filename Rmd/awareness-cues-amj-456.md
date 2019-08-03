---
title: "awareness-cues-amj-456"
author: "Stephen Downing"
date: "7/3/2019"
output: 
  html_document: 
    keep_md: yes
---



## R Markdown

## Awareness Cues AMJ Parts 4,5,6

Markdown for Parts 4,5,6 of Awareness Cues computations. 


```r
library(btergm)
```

```
## Loading required package: xergm.common
```

```
## Loading required package: ergm
```

```
## Loading required package: network
```

```
## network: Classes for Relational Data
## Version 1.15 created on 2019-04-01.
## copyright (c) 2005, Carter T. Butts, University of California-Irvine
##                     Mark S. Handcock, University of California -- Los Angeles
##                     David R. Hunter, Penn State University
##                     Martina Morris, University of Washington
##                     Skye Bender-deMoll, University of Washington
##  For citation information, type citation("network").
##  Type help("network-package") to get started.
```

```
## 
## ergm: version 3.10.1, created on 2019-05-14
## Copyright (c) 2019, Mark S. Handcock, University of California -- Los Angeles
##                     David R. Hunter, Penn State University
##                     Carter T. Butts, University of California -- Irvine
##                     Steven M. Goodreau, University of Washington
##                     Pavel N. Krivitsky, University of Wollongong
##                     Martina Morris, University of Washington
##                     with contributions from
##                     Li Wang
##                     Kirk Li, University of Washington
##                     Skye Bender-deMoll, University of Washington
##                     Chad Klumb
## Based on "statnet" project software (statnet.org).
## For license and citation information see statnet.org/attribution
## or type citation("ergm").
```

```
## NOTE: Versions before 3.6.1 had a bug in the implementation of the
## bd() constriant which distorted the sampled distribution somewhat.
## In addition, Sampson's Monks datasets had mislabeled vertices. See
## the NEWS and the documentation for more details.
```

```
## NOTE: Some common term arguments pertaining to vertex attribute
## and level selection have changed in 3.10.0. See terms help for
## more details. Use 'options(ergm.term=list(version="3.9.4"))' to
## use old behavior.
```

```
## 
## Attaching package: 'xergm.common'
```

```
## The following object is masked from 'package:ergm':
## 
##     gof
```

```
## Loading required package: ggplot2
```

```
## Registered S3 methods overwritten by 'ggplot2':
##   method         from 
##   [.quosures     rlang
##   c.quosures     rlang
##   print.quosures rlang
```

```
## Registered S3 methods overwritten by 'btergm':
##   method    from
##   print.gof ergm
##   plot.gof  ergm
```

```
## Package:  btergm
## Version:  1.9.4
## Date:     2019-05-12
## Authors:  Philip Leifeld (University of Essex)
##           Skyler J. Cranmer (The Ohio State University)
##           Bruce A. Desmarais (Pennsylvania State University)
```

```r
library(parallel)
library(texreg)
```

```
## Version:  1.36.24
## Date:     2019-04-07
## Author:   Philip Leifeld (University of Essex)
## 
## Please cite the JSS article in your publications -- see citation("texreg").
```

```r
## DIRECTORIES
data_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/awareness-cues-amj"
work_dir <- data_dir
```




```r
##===============================
## set analysis params
##-------------------------------
firm_i <- 'qualtrics'
d <- 2
ncpus <- 4
parallel <- "multicore"
nPeriods <- 11  ## 5
##-------------------------------

data_file <- file.path(data_dir,sprintf('%s_d%s.rds',firm_i,d))
nets <- readRDS(data_file)

#cache lists
if (nPeriods < length(nets))   
  nets <- nets[(length(nets)-nPeriods+1):length(nets)] 
```



```r
cat("\n------------ estimating TERGM for:",firm_i,'--------------\n')
```

```
## 
## ------------ estimating TERGM for: qualtrics --------------
```

```r
cat(sprintf("Using %s cores\n", detectCores()))
```

```
## Using 4 cores
```

```r
## make MMC nets list
mmc <- lapply(nets, function(net) as.matrix(net %n% 'mmc'))
cpc <- lapply(nets, function(net) as.matrix(net %n% 'coop'))
cpp <- lapply(nets, function(net) as.matrix(net %n% 'coop_past'))
cpa <- lapply(nets, function(net) as.matrix(net %n% 'coop') + as.matrix(net %n% 'coop_past') )
cossim <- lapply(nets, function(net) as.matrix(net %n% 'cat_cos_sim'))
centjoin <- lapply(nets, function(net) as.matrix(net %n% 'joint_cent_pow_n0_4'))
centratio <- lapply(nets, function(net) as.matrix(net %n% 'cent_ratio_pow_n0_4'))
shcomp <- lapply(nets, function(net) as.matrix(net %n% 'shared_competitor')) 
shinv <- lapply(nets, function(net) as.matrix(net %n% 'shared_investor_nd'))

####################### DEFINE MODELS ###################################

m4 <-   nets ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed=T) + 
  nodematch("ipo_status", diff = F) + 
  nodematch("state_code", diff = F) + 
  nodecov("age") + absdiff("age") + 
  nodecov("employee_na_age") +
  nodecov("sales_na_0_mn") +
  edgecov(cossim) +
  edgecov(shinv) +   
  edgecov(mmc) + 
  memory(type = "stability", lag = 1) + 
  timecov(transform = function(t) t) +
  nodecov("genidx_multilevel") + 
  nodecov("cent_pow_n0_4") + absdiff("cent_pow_n0_4") + 
  cycle(3) + cycle(4) + cycle(5) 

################################ end models#######################

##
# DEFINE MODEL and MODEL NAME TO COMPUTE
## 
m_x <- 'm4'
##

##==============================================
## BTERGM
##----------------------------------------------

# SET RESAMPLES
##
R <- 100

## RUN TERGM
fit4 <- btergm(get(m_x), R=R, parallel = parallel, ncpus = ncpus)
```

```
## Mean transformed timecov values:
```

```
##  t=1 t=2 t=3 t=4 t=5 t=6 t=7 t=8 t=9 t=10 t=11
##    1   2   3   4   5   6   7   8   9   10   11
```

```
## 
## Initial dimensions of the network and covariates:
```

```
##                t=2 t=3 t=4 t=5 t=6 t=7 t=8 t=9 t=10 t=11
## nets (row)     168 168 168 168 168 168 168 168  168  168
## nets (col)     168 168 168 168 168 168 168 168  168  168
## cossim (row)   168 168 168 168 168 168 168 168  168  168
## cossim (col)   168 168 168 168 168 168 168 168  168  168
## shinv (row)    168 168 168 168 168 168 168 168  168  168
## shinv (col)    168 168 168 168 168 168 168 168  168  168
## mmc (row)      168 168 168 168 168 168 168 168  168  168
## mmc (col)      168 168 168 168 168 168 168 168  168  168
## memory (row)   168 168 168 168 168 168 168 168  168  168
## memory (col)   168 168 168 168 168 168 168 168  168  168
## timecov1 (row) 168 168 168 168 168 168 168 168  168  168
## timecov1 (col) 168 168 168 168 168 168 168 168  168  168
```

```
## 
## All networks are conformable.
```

```
## 
## Dimensions of the network and covariates after adjustment:
```

```
##                t=2 t=3 t=4 t=5 t=6 t=7 t=8 t=9 t=10 t=11
## nets (row)     168 168 168 168 168 168 168 168  168  168
## nets (col)     168 168 168 168 168 168 168 168  168  168
## cossim (row)   168 168 168 168 168 168 168 168  168  168
## cossim (col)   168 168 168 168 168 168 168 168  168  168
## shinv (row)    168 168 168 168 168 168 168 168  168  168
## shinv (col)    168 168 168 168 168 168 168 168  168  168
## mmc (row)      168 168 168 168 168 168 168 168  168  168
## mmc (col)      168 168 168 168 168 168 168 168  168  168
## memory (row)   168 168 168 168 168 168 168 168  168  168
## memory (col)   168 168 168 168 168 168 168 168  168  168
## timecov1 (row) 168 168 168 168 168 168 168 168  168  168
## timecov1 (col) 168 168 168 168 168 168 168 168  168  168
```

```
## 
## Starting pseudolikelihood estimation with 100 bootstrapping replications using multicore forking on 4 cores...
```

```
## Done.
```

```r
## SAVE SERIALIZED
fits.file <- sprintf('fit_%s_pd%s_R%s_%s_d%s.rds', firm_i, nPeriods, R, m_x, d)
saveRDS(fit4, file=file.path(data_dir,fits.file))
## SAVE FORMATTED REGRESSION TABLE
html.file <- sprintf('%s_tergm_results_pd%s_R%s_%s_d%s.html',  firm_i, nPeriods, R, m_x, d)
htmlreg(fit4, digits = 2, file=file.path(data_dir,html.file))
```

```
## The table was written to the file 'C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/awareness-cues-amj/qualtrics_tergm_results_pd11_R100_m4_d2.html'.
```



```r
##================================================
## Goodness of Fit
##------------------------------------------------
gof4 <- gof(fit4, statistics=c(dsp, esp, deg, geodesic), nsim=10)
```

```
## 
## Starting GOF assessment on a single computing core....
```

```
## Mean transformed timecov values:
```

```
##  t=1 t=2 t=3 t=4 t=5 t=6 t=7 t=8 t=9 t=10 t=11
##    1   2   3   4   5   6   7   8   9   10   11
```

```
## 
## Initial dimensions of the network and covariates:
```

```
##                t=2 t=3 t=4 t=5 t=6 t=7 t=8 t=9 t=10 t=11
## nets (row)     168 168 168 168 168 168 168 168  168  168
## nets (col)     168 168 168 168 168 168 168 168  168  168
## cossim (row)   168 168 168 168 168 168 168 168  168  168
## cossim (col)   168 168 168 168 168 168 168 168  168  168
## shinv (row)    168 168 168 168 168 168 168 168  168  168
## shinv (col)    168 168 168 168 168 168 168 168  168  168
## mmc (row)      168 168 168 168 168 168 168 168  168  168
## mmc (col)      168 168 168 168 168 168 168 168  168  168
## memory (row)   168 168 168 168 168 168 168 168  168  168
## memory (col)   168 168 168 168 168 168 168 168  168  168
## timecov1 (row) 168 168 168 168 168 168 168 168  168  168
## timecov1 (col) 168 168 168 168 168 168 168 168  168  168
```

```
## 
## All networks are conformable.
```

```
## 
## Dimensions of the network and covariates after adjustment:
```

```
##                t=2 t=3 t=4 t=5 t=6 t=7 t=8 t=9 t=10 t=11
## nets (row)     168 168 168 168 168 168 168 168  168  168
## nets (col)     168 168 168 168 168 168 168 168  168  168
## cossim (row)   168 168 168 168 168 168 168 168  168  168
## cossim (col)   168 168 168 168 168 168 168 168  168  168
## shinv (row)    168 168 168 168 168 168 168 168  168  168
## shinv (col)    168 168 168 168 168 168 168 168  168  168
## mmc (row)      168 168 168 168 168 168 168 168  168  168
## mmc (col)      168 168 168 168 168 168 168 168  168  168
## memory (row)   168 168 168 168 168 168 168 168  168  168
## memory (col)   168 168 168 168 168 168 168 168  168  168
## timecov1 (row) 168 168 168 168 168 168 168 168  168  168
## timecov1 (col) 168 168 168 168 168 168 168 168  168  168
```

```
## 
## No 'target' network(s) provided. Using networks on the left-hand side of the model formula as observed networks.
```

```
## Simulating 10 networks from the following formula:
##  nets[[1]] ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + nodematch("ipo_status", diff = F) + nodematch("state_code", diff = F) + nodecov("age") + absdiff("age") + nodecov("employee_na_age") + nodecov("sales_na_0_mn") + edgecov(cossim[[1]]) + edgecov(shinv[[1]]) + edgecov(mmc[[1]]) + edgecov(memory[[1]]) + edgecov(timecov1[[1]]) + nodecov("genidx_multilevel") + nodecov("cent_pow_n0_4") + absdiff("cent_pow_n0_4") + cycle(3) + cycle(4) + cycle(5)
```

```
## Simulating 10 networks from the following formula:
##  nets[[2]] ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + nodematch("ipo_status", diff = F) + nodematch("state_code", diff = F) + nodecov("age") + absdiff("age") + nodecov("employee_na_age") + nodecov("sales_na_0_mn") + edgecov(cossim[[2]]) + edgecov(shinv[[2]]) + edgecov(mmc[[2]]) + edgecov(memory[[2]]) + edgecov(timecov1[[2]]) + nodecov("genidx_multilevel") + nodecov("cent_pow_n0_4") + absdiff("cent_pow_n0_4") + cycle(3) + cycle(4) + cycle(5)
```

```
## Simulating 10 networks from the following formula:
##  nets[[3]] ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + nodematch("ipo_status", diff = F) + nodematch("state_code", diff = F) + nodecov("age") + absdiff("age") + nodecov("employee_na_age") + nodecov("sales_na_0_mn") + edgecov(cossim[[3]]) + edgecov(shinv[[3]]) + edgecov(mmc[[3]]) + edgecov(memory[[3]]) + edgecov(timecov1[[3]]) + nodecov("genidx_multilevel") + nodecov("cent_pow_n0_4") + absdiff("cent_pow_n0_4") + cycle(3) + cycle(4) + cycle(5)
```

```
## Simulating 10 networks from the following formula:
##  nets[[4]] ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + nodematch("ipo_status", diff = F) + nodematch("state_code", diff = F) + nodecov("age") + absdiff("age") + nodecov("employee_na_age") + nodecov("sales_na_0_mn") + edgecov(cossim[[4]]) + edgecov(shinv[[4]]) + edgecov(mmc[[4]]) + edgecov(memory[[4]]) + edgecov(timecov1[[4]]) + nodecov("genidx_multilevel") + nodecov("cent_pow_n0_4") + absdiff("cent_pow_n0_4") + cycle(3) + cycle(4) + cycle(5)
```

```
## Simulating 10 networks from the following formula:
##  nets[[5]] ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + nodematch("ipo_status", diff = F) + nodematch("state_code", diff = F) + nodecov("age") + absdiff("age") + nodecov("employee_na_age") + nodecov("sales_na_0_mn") + edgecov(cossim[[5]]) + edgecov(shinv[[5]]) + edgecov(mmc[[5]]) + edgecov(memory[[5]]) + edgecov(timecov1[[5]]) + nodecov("genidx_multilevel") + nodecov("cent_pow_n0_4") + absdiff("cent_pow_n0_4") + cycle(3) + cycle(4) + cycle(5)
```

```
## Simulating 10 networks from the following formula:
##  nets[[6]] ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + nodematch("ipo_status", diff = F) + nodematch("state_code", diff = F) + nodecov("age") + absdiff("age") + nodecov("employee_na_age") + nodecov("sales_na_0_mn") + edgecov(cossim[[6]]) + edgecov(shinv[[6]]) + edgecov(mmc[[6]]) + edgecov(memory[[6]]) + edgecov(timecov1[[6]]) + nodecov("genidx_multilevel") + nodecov("cent_pow_n0_4") + absdiff("cent_pow_n0_4") + cycle(3) + cycle(4) + cycle(5)
```

```
## Simulating 10 networks from the following formula:
##  nets[[7]] ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + nodematch("ipo_status", diff = F) + nodematch("state_code", diff = F) + nodecov("age") + absdiff("age") + nodecov("employee_na_age") + nodecov("sales_na_0_mn") + edgecov(cossim[[7]]) + edgecov(shinv[[7]]) + edgecov(mmc[[7]]) + edgecov(memory[[7]]) + edgecov(timecov1[[7]]) + nodecov("genidx_multilevel") + nodecov("cent_pow_n0_4") + absdiff("cent_pow_n0_4") + cycle(3) + cycle(4) + cycle(5)
```

```
## Simulating 10 networks from the following formula:
##  nets[[8]] ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + nodematch("ipo_status", diff = F) + nodematch("state_code", diff = F) + nodecov("age") + absdiff("age") + nodecov("employee_na_age") + nodecov("sales_na_0_mn") + edgecov(cossim[[8]]) + edgecov(shinv[[8]]) + edgecov(mmc[[8]]) + edgecov(memory[[8]]) + edgecov(timecov1[[8]]) + nodecov("genidx_multilevel") + nodecov("cent_pow_n0_4") + absdiff("cent_pow_n0_4") + cycle(3) + cycle(4) + cycle(5)
```

```
## Simulating 10 networks from the following formula:
##  nets[[9]] ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + nodematch("ipo_status", diff = F) + nodematch("state_code", diff = F) + nodecov("age") + absdiff("age") + nodecov("employee_na_age") + nodecov("sales_na_0_mn") + edgecov(cossim[[9]]) + edgecov(shinv[[9]]) + edgecov(mmc[[9]]) + edgecov(memory[[9]]) + edgecov(timecov1[[9]]) + nodecov("genidx_multilevel") + nodecov("cent_pow_n0_4") + absdiff("cent_pow_n0_4") + cycle(3) + cycle(4) + cycle(5)
```

```
## Simulating 10 networks from the following formula:
##  nets[[10]] ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + nodematch("ipo_status", diff = F) + nodematch("state_code", diff = F) + nodecov("age") + absdiff("age") + nodecov("employee_na_age") + nodecov("sales_na_0_mn") + edgecov(cossim[[10]]) + edgecov(shinv[[10]]) + edgecov(mmc[[10]]) + edgecov(memory[[10]]) + edgecov(timecov1[[10]]) + nodecov("genidx_multilevel") + nodecov("cent_pow_n0_4") + absdiff("cent_pow_n0_4") + cycle(3) + cycle(4) + cycle(5)
```

```
## 10 networks from which simulations are drawn were provided.
```

```
## Processing statistic: Dyad-wise shared partners
```

```
## Processing statistic: Edge-wise shared partners
```

```
## Processing statistic: Degree
```

```
## Processing statistic: Geodesic distances
```

```r
print(gof4)
```

```
## Dyad-wise shared partners
```

```
##    obs: mean median   min   max sim: mean median   min   max    Pr(>z)    
## 0    24614.4  25452 20730 27614  24594.54  25885 19942 27512 0.9840922    
## 1     2568.4   1947   288  5558   2511.26   1401   380  6042 0.9401320    
## 2      618.4    460   118  1320    688.66    573   106  1710 0.6788767    
## 3      147.4    111    30   290    151.02    115    34   348 0.9229278    
## 4       63.2     59     4   130     65.26     49     6   154 0.9044678    
## 5       24.6     16     0    66     24.44     14     0    66 0.9856136    
## 6        8.0      4     0    24      8.70      4     0    32 0.8302252    
## 7        2.8      2     0     8      3.94      2     0    14 0.3156949    
## 8        2.4      1     0     8      1.80      0     0     8 0.5639004    
## 9        2.2      0     0     6      1.76      0     0    10 0.6507504    
## 10       0.8      0     0     2      1.04      0     0     8 0.5244368    
## 11       0.6      0     0     2      0.58      0     0     4 0.9513169    
## 12       0.2      0     0     2      0.34      0     0     4 0.5283204    
## 13       0.0      0     0     0      0.28      0     0     4 0.0003334 ***
## 14       0.6      0     0     2      0.24      0     0     2 0.2764115    
## 15       0.0      0     0     0      0.36      0     0     4 0.0001570 ***
## 16       0.0      0     0     0      0.06      0     0     2 0.0832486 .  
## 17       0.0      0     0     0      0.24      0     0     2 0.0003873 ***
## 18       0.4      0     0     2      0.06      0     0     2 0.2367899    
## 19       0.0      0     0     0      0.02      0     0     2 0.3197485    
## 20       0.0      0     0     0      0.00      0     0     0 1.0000000    
## 21       0.4      0     0     2      0.00      0     0     0 0.1678507    
## 22       0.4      0     0     2      0.00      0     0     0 0.1678507    
## 23       0.6      0     0     2      0.22      0     0     2 0.2516930    
## 24       0.2      0     0     2      0.40      0     0     2 0.3716113    
## 25       0.0      0     0     0      0.20      0     0     2 0.0012748 ** 
## 26       0.0      0     0     0      0.20      0     0     2 0.0012748 ** 
## 27       0.0      0     0     0      0.14      0     0     2 0.0075026 ** 
## 28       0.0      0     0     0      0.12      0     0     2 0.0135568 *  
## 29       0.0      0     0     0      0.06      0     0     2 0.0832486 .  
## 30       0.0      0     0     0      0.06      0     0     2 0.0832486 .  
## 31       0.0      0     0     0      0.00      0     0     0 1.0000000
```

```
## 
## Note: Small p-values indicate a significant difference 
##       between simulations and observed network(s).
```

```
## Edge-wise shared partners
```

```
##    obs: mean median min max sim: mean median min max    Pr(>z)    
## 0      102.2    107  28 168     96.40     97  34 160 0.7340861    
## 1      130.6    112  62 218    132.88    117  46 266 0.9128038    
## 2       65.2     60  22 124     63.36     56  14 126 0.8903381    
## 3       51.2     43  16  90     47.30     41  18  94 0.6859690    
## 4       24.4     24   0  58     29.74     25   2  80 0.4626491    
## 5       12.6      7   0  32     12.06      8   0  34 0.9053355    
## 6        4.6      3   0  12      5.82      2   0  18 0.4899202    
## 7        2.4      1   0   8      3.10      2   0  10 0.5140188    
## 8        1.8      1   0   6      1.26      0   0   6 0.4676324    
## 9        2.2      0   0   6      1.62      0   0  10 0.5516448    
## 10       0.8      0   0   2      0.86      0   0   8 0.8701063    
## 11       0.4      0   0   2      0.54      0   0   2 0.6282835    
## 12       0.0      0   0   0      0.26      0   0   2 0.0002124 ***
## 13       0.0      0   0   0      0.14      0   0   2 0.0075026 ** 
## 14       0.0      0   0   0      0.04      0   0   2 0.1583399    
## 15       0.0      0   0   0      0.20      0   0   4 0.0068956 ** 
## 16       0.0      0   0   0      0.04      0   0   2 0.1583399    
## 17       0.0      0   0   0      0.20      0   0   2 0.0012748 ** 
## 18       0.4      0   0   2      0.06      0   0   2 0.2367899    
## 19       0.0      0   0   0      0.02      0   0   2 0.3197485    
## 20       0.0      0   0   0      0.00      0   0   0 1.0000000    
## 21       0.4      0   0   2      0.00      0   0   0 0.1678507    
## 22       0.4      0   0   2      0.00      0   0   0 0.1678507    
## 23       0.6      0   0   2      0.22      0   0   2 0.2516930    
## 24       0.2      0   0   2      0.40      0   0   2 0.3716113    
## 25       0.0      0   0   0      0.20      0   0   2 0.0012748 ** 
## 26       0.0      0   0   0      0.20      0   0   2 0.0012748 ** 
## 27       0.0      0   0   0      0.14      0   0   2 0.0075026 ** 
## 28       0.0      0   0   0      0.12      0   0   2 0.0135568 *  
## 29       0.0      0   0   0      0.06      0   0   2 0.0832486 .  
## 30       0.0      0   0   0      0.06      0   0   2 0.0832486 .  
## 31       0.0      0   0   0      0.00      0   0   0 1.0000000
```

```
## 
## Note: Small p-values indicate a significant difference 
##       between simulations and observed network(s).
```

```
## Degree
```

```
##    obs: mean median min max sim: mean median min max    Pr(>z)    
## 0       73.1   79.5  16 133     76.15   80.5  14 128 0.8431425    
## 1       30.2   28.5   7  53     27.65   26.0  11  49 0.6572661    
## 2       17.1   15.5  10  25     17.75   15.0   8  37 0.7675791    
## 3       12.4   12.5   4  21     11.92   11.0   1  28 0.8319513    
## 4        7.8    7.5   2  14      7.67    8.0   1  19 0.9320433    
## 5        8.4    8.0   5  13      7.25    6.0   1  16 0.2650948    
## 6        6.1    5.0   2  11      5.51    4.0   0  14 0.5882696    
## 7        2.2    3.0   0   4      2.99    3.0   0   7 0.1215835    
## 8        2.7    3.0   1   4      2.48    2.0   0   8 0.5955534    
## 9        1.1    1.0   0   3      1.53    1.5   0   5 0.1787751    
## 10       0.9    1.0   0   2      1.10    1.0   0   3 0.5095212    
## 11       0.9    1.0   0   2      0.58    0.0   0   2 0.1227472    
## 12       1.5    1.0   0   3      0.96    1.0   0   5 0.1623628    
## 13       0.1    0.0   0   1      0.65    1.0   0   2 0.0002273 ***
## 14       0.1    0.0   0   1      0.52    0.0   0   3 0.0024280 ** 
## 15       0.4    0.0   0   1      0.38    0.0   0   2 0.9098874    
## 16       0.0    0.0   0   0      0.22    0.0   0   2  6.65e-06 ***
## 17       0.1    0.0   0   1      0.18    0.0   0   1 0.4700225    
## 18       0.1    0.0   0   1      0.12    0.0   0   1 0.8526768    
## 19       0.2    0.0   0   1      0.13    0.0   0   2 0.6232612    
## 20       0.3    0.0   0   1      0.17    0.0   0   1 0.4277267    
## 21       0.3    0.0   0   1      0.18    0.0   0   2 0.4665330    
## 22       0.0    0.0   0   0      0.10    0.0   0   2 0.0068956 ** 
## 23       0.0    0.0   0   0      0.06    0.0   0   1 0.0135568 *  
## 24       0.1    0.0   0   1      0.03    0.0   0   1 0.5066781    
## 25       0.1    0.0   0   1      0.11    0.0   0   1 0.9257392    
## 26       0.2    0.0   0   1      0.17    0.0   0   2 0.8335314    
## 27       0.2    0.0   0   1      0.07    0.0   0   1 0.3616470    
## 28       0.0    0.0   0   0      0.08    0.0   0   1 0.0041580 ** 
## 29       0.1    0.0   0   1      0.01    0.0   0   1 0.3933760    
## 30       0.0    0.0   0   0      0.02    0.0   0   1 0.1583399    
## 31       0.0    0.0   0   0      0.02    0.0   0   1 0.1583399    
## 32       0.0    0.0   0   0      0.05    0.0   0   1 0.0245895 *  
## 33       0.0    0.0   0   0      0.03    0.0   0   1 0.0832486 .  
## 34       0.0    0.0   0   0      0.01    0.0   0   1 0.3197485    
## 35       0.0    0.0   0   0      0.03    0.0   0   1 0.0832486 .  
## 36       0.0    0.0   0   0      0.02    0.0   0   1 0.1583399    
## 37       0.1    0.0   0   1      0.03    0.0   0   1 0.5066781    
## 38       0.0    0.0   0   0      0.01    0.0   0   1 0.3197485    
## 39       0.0    0.0   0   0      0.05    0.0   0   1 0.0245895 *  
## 40       0.0    0.0   0   0      0.03    0.0   0   1 0.0832486 .  
## 41       0.1    0.0   0   1      0.03    0.0   0   1 0.5066781    
## 42       0.1    0.0   0   1      0.03    0.0   0   1 0.5066781    
## 43       0.1    0.0   0   1      0.03    0.0   0   1 0.5066781    
## 44       0.2    0.0   0   1      0.01    0.0   0   1 0.1886567    
## 45       0.2    0.0   0   1      0.05    0.0   0   1 0.2942787    
## 46       0.1    0.0   0   1      0.05    0.0   0   1 0.6359084    
## 47       0.1    0.0   0   1      0.06    0.0   0   1 0.7053387    
## 48       0.1    0.0   0   1      0.06    0.0   0   1 0.7053387    
## 49       0.0    0.0   0   0      0.11    0.0   0   1 0.0007038 ***
## 50       0.0    0.0   0   0      0.12    0.0   0   2 0.0023040 ** 
## 51       0.1    0.0   0   1      0.15    0.0   0   2 0.6543701    
## 52       0.1    0.0   0   1      0.07    0.0   0   1 0.7771715    
## 53       0.0    0.0   0   0      0.08    0.0   0   1 0.0041580 ** 
## 54       0.0    0.0   0   0      0.05    0.0   0   1 0.0245895 *  
## 55       0.0    0.0   0   0      0.01    0.0   0   1 0.3197485    
## 56       0.0    0.0   0   0      0.01    0.0   0   1 0.3197485    
## 57       0.0    0.0   0   0      0.02    0.0   0   1 0.1583399    
## 58       0.0    0.0   0   0      0.01    0.0   0   1 0.3197485    
## 59       0.0    0.0   0   0      0.03    0.0   0   1 0.0832486 .  
## 60       0.0    0.0   0   0      0.00    0.0   0   0 1.0000000
```

```
## 
## Note: Small p-values indicate a significant difference 
##       between simulations and observed network(s).
```

```
## Geodesic distances
```

```
##     obs: mean median  min   max sim: mean median  min   max  Pr(>z)  
## 1       400.4    359  132   682    397.30    365  142   718 0.96725  
## 2      3143.4   2352  338  6812   3160.56   1905  436  7560 0.98539  
## 3      3974.8   2213  388  9568   3733.28   2434  524  9910 0.85351  
## 4      2770.6   2097  256  5824   2038.74   1665  168  4990 0.33347  
## 5       319.8    271   68   660    326.40    342   20   926 0.92195  
## 6        52.6     26    0   196     58.12     23    0   348 0.81436  
## 7        10.0      0    0    36      9.08      0    0   126 0.86934  
## 8         0.0      0    0     0      0.66      0    0    24 0.02432 *
## 9         0.0      0    0     0      0.00      0    0     0 1.00000  
## Inf   17384.4  20450 5104 26866  18331.86  20657 4800 26642 0.75319
```

```
## 
## Note: Small p-values indicate a significant difference 
##       between simulations and observed network(s).
```

```r
## plot
plot(gof4)
```

![](awareness-cues-amj-456_files/figure-html/gof-1.png)<!-- -->




```r
##================================================
## COMPARE ESTIMATION ALGORITHM
##------------------------------------------------

set.seed(1111)
## estimate the TERGM with bootstrapped PMLE
fit4mc <- mtergm(get(m_x), ctrl=control.ergm(seed = 1111))
```

```
## Mean transformed timecov values:
```

```
##  t=1 t=2 t=3 t=4 t=5 t=6 t=7 t=8 t=9 t=10 t=11
##    1   2   3   4   5   6   7   8   9   10   11
```

```
## 
## Initial dimensions of the network and covariates:
```

```
##                t=2 t=3 t=4 t=5 t=6 t=7 t=8 t=9 t=10 t=11
## nets (row)     168 168 168 168 168 168 168 168  168  168
## nets (col)     168 168 168 168 168 168 168 168  168  168
## cossim (row)   168 168 168 168 168 168 168 168  168  168
## cossim (col)   168 168 168 168 168 168 168 168  168  168
## shinv (row)    168 168 168 168 168 168 168 168  168  168
## shinv (col)    168 168 168 168 168 168 168 168  168  168
## mmc (row)      168 168 168 168 168 168 168 168  168  168
## mmc (col)      168 168 168 168 168 168 168 168  168  168
## memory (row)   168 168 168 168 168 168 168 168  168  168
## memory (col)   168 168 168 168 168 168 168 168  168  168
## timecov1 (row) 168 168 168 168 168 168 168 168  168  168
## timecov1 (col) 168 168 168 168 168 168 168 168  168  168
```

```
## 
## All networks are conformable.
```

```
## 
## Dimensions of the network and covariates after adjustment:
```

```
##                t=2 t=3 t=4 t=5 t=6 t=7 t=8 t=9 t=10 t=11
## nets (row)     168 168 168 168 168 168 168 168  168  168
## nets (col)     168 168 168 168 168 168 168 168  168  168
## cossim (row)   168 168 168 168 168 168 168 168  168  168
## cossim (col)   168 168 168 168 168 168 168 168  168  168
## shinv (row)    168 168 168 168 168 168 168 168  168  168
## shinv (col)    168 168 168 168 168 168 168 168  168  168
## mmc (row)      168 168 168 168 168 168 168 168  168  168
## mmc (col)      168 168 168 168 168 168 168 168  168  168
## memory (row)   168 168 168 168 168 168 168 168  168  168
## memory (col)   168 168 168 168 168 168 168 168  168  168
## timecov1 (row) 168 168 168 168 168 168 168 168  168  168
## timecov1 (col) 168 168 168 168 168 168 168 168  168  168
```

```
## Estimating...
```

```
## Starting maximum pseudolikelihood estimation (MPLE):
```

```
## Evaluating the predictor and response matrix.
```

```
## Maximizing the pseudolikelihood.
```

```
## Finished MPLE.
```

```
## Starting Monte Carlo maximum likelihood estimation (MCMLE):
```

```
## Iteration 1 of at most 20:
```

```
## Optimizing with step length 0.263988225860428.
```

```
## The log-likelihood improved by 2.846.
```

```
## Iteration 2 of at most 20:
```

```
## Optimizing with step length 0.118130424516662.
```

```
## The log-likelihood improved by 1.997.
```

```
## Iteration 3 of at most 20:
```

```
## Optimizing with step length 0.137912842628196.
```

```
## The log-likelihood improved by 1.867.
```

```
## Iteration 4 of at most 20:
```

```
## Optimizing with step length 0.185868846371085.
```

```
## The log-likelihood improved by 1.66.
```

```
## Iteration 5 of at most 20:
```

```
## Optimizing with step length 0.0618567092390864.
```

```
## The log-likelihood improved by 2.722.
```

```
## Iteration 6 of at most 20:
```

```
## Optimizing with step length 0.183692905265568.
```

```
## The log-likelihood improved by 2.334.
```

```
## Iteration 7 of at most 20:
```

```
## Optimizing with step length 0.0521149937336016.
```

```
## The log-likelihood improved by 2.059.
```

```
## Iteration 8 of at most 20:
```

```
## Optimizing with step length 0.0219064701676662.
```

```
## The log-likelihood improved by 1.452.
```

```
## Iteration 9 of at most 20:
```

```
## Optimizing with step length 0.153135590802792.
```

```
## The log-likelihood improved by 2.324.
```

```
## Iteration 10 of at most 20:
```

```
## Optimizing with step length 0.137134837199849.
```

```
## The log-likelihood improved by 1.445.
```

```
## Iteration 11 of at most 20:
```

```
## Optimizing with step length 0.247173407130022.
```

```
## The log-likelihood improved by 1.711.
```

```
## Iteration 12 of at most 20:
```

```
## Optimizing with step length 0.256508177175332.
```

```
## The log-likelihood improved by 1.562.
```

```
## Iteration 13 of at most 20:
```

```
## Optimizing with step length 0.0657523307296441.
```

```
## The log-likelihood improved by 1.474.
```

```
## Iteration 14 of at most 20:
```

```
## Optimizing with step length 0.13890694972477.
```

```
## The log-likelihood improved by 2.173.
```

```
## Iteration 15 of at most 20:
```

```
## Optimizing with step length 0.181754127787292.
```

```
## The log-likelihood improved by 1.923.
```

```
## Iteration 16 of at most 20:
```

```
## Optimizing with step length 0.0471135041644797.
```

```
## The log-likelihood improved by 1.636.
```

```
## Iteration 17 of at most 20:
```

```
## Optimizing with step length 0.0869250395344613.
```

```
## The log-likelihood improved by 1.785.
```

```
## Iteration 18 of at most 20:
```

```
## Optimizing with step length 0.160340839388494.
```

```
## The log-likelihood improved by 1.676.
```

```
## Iteration 19 of at most 20:
```

```
## Optimizing with step length 0.355943293387886.
```

```
## The log-likelihood improved by 2.003.
```

```
## Iteration 20 of at most 20:
```

```
## Optimizing with step length 0.10322767499567.
```

```
## The log-likelihood improved by 1.813.
```

```
## MCMLE estimation did not converge after 20 iterations. The estimated coefficients may not be accurate. Estimation may be resumed by passing the coefficients as initial values; see 'init' under ?control.ergm for details.
```

```
## Finished MCMLE.
```

```
## Evaluating log-likelihood at the estimate. Using 20 bridges: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 .
## This model was fit using MCMC.  To examine model diagnostics and
## check for degeneracy, use the mcmc.diagnostics() function.
## Done.
```

```r
## SAVE SERIALIZED
fits.file <- sprintf('fit_%s_pd%s_R%s_%s_d%s.rds', firm_i, nPeriods, R, m_x, d)
saveRDS(fit4mc, file=file.path(data_dir,fits.file))
## SAVE FORMATTED REGRESSION TABLE
html.file <- sprintf('%s_tergm_results_pd%s_R%s_%s_d%s.html',  firm_i, nPeriods, R, m_x, d)
htmlreg(fit4mc, digits = 2, file=file.path(data_dir,html.file))
```

```
## The table was written to the file 'C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/awareness-cues-amj/qualtrics_tergm_results_pd11_R100_m4_d2.html'.
```

```r
cat('finished successfully.')
```

```
## finished successfully.
```

```r
## COMPARISON TABLE 
## Cache model fits list
fits <- list(PMLE=fit4, MCMCMLE=fit4mc)
## Echo model comparison table to screen
screenreg(fits, digits = 3)
```

```
## 
## ===========================================================
##                            PMLE              MCMCMLE       
## -----------------------------------------------------------
## edges                          -3.251 *          -4.459 ***
##                            [-5.654; -0.935]      (0.262)   
## gwesp.fixed.0                   0.358 *           0.424 ** 
##                            [ 0.057;  0.572]      (0.161)   
## gwdegree                       -0.907                      
##                            [-1.950;  0.074]                
## nodematch.ipo_status            0.243             0.691    
##                            [-1.033;  2.623]      (0.440)   
## nodematch.state_code            0.025             0.006    
##                            [-0.339;  0.539]      (0.220)   
## nodecov.age                    -0.105 *          -0.106 ***
##                            [-0.164; -0.065]      (0.017)   
## absdiff.age                     0.124 *           0.124 ***
##                            [ 0.083;  0.186]      (0.018)   
## nodecov.employee_na_age         0.000             0.000    
##                            [-0.000;  0.001]      (0.000)   
## nodecov.sales_na_0_mn          -0.000            -0.000    
##                            [-0.005;  0.003]      (0.001)   
## edgecov.cossim[[i]]             1.835 *                    
##                            [ 1.335;  2.135]                
## edgecov.shinv[[i]]              0.068                      
##                            [-0.106;  0.256]                
## edgecov.mmc[[i]]               -3.051                      
##                            [-5.760;  0.804]                
## edgecov.memory[[i]]             4.925 *                    
##                            [ 4.691;  6.095]                
## edgecov.timecov1[[i]]           0.033                      
##                            [-0.187;  0.316]                
## nodecov.genidx_multilevel       1.776 *           2.156 ***
##                            [ 1.077;  2.220]      (0.125)   
## nodecov.cent_pow_n0_4          -0.034            -0.009    
##                            [-0.111;  0.045]      (0.039)   
## absdiff.cent_pow_n0_4           0.272 *           0.323 ***
##                            [ 0.143;  0.455]      (0.061)   
## cycle3                          0.483 *           0.321    
##                            [ 0.235;  0.702]      (0.200)   
## cycle4                          0.123 *           0.053    
##                            [ 0.086;  0.231]      (0.027)   
## cycle5                         -0.026 *          -0.021 ***
##                            [-0.053; -0.019]      (0.004)   
## gwdeg.fixed.0                                    -0.423 *  
##                                                  (0.204)   
## edgecov.cossim                                    2.524 ***
##                                                  (0.296)   
## edgecov.shinv                                    -0.272    
##                                                  (0.192)   
## edgecov.mmc                                      -0.246    
##                                                  (0.235)   
## edgecov.memory                                    4.840 ***
##                                                  (0.521)   
## edgecov.timecov1                                  0.080 *  
##                                                  (0.033)   
## -----------------------------------------------------------
## Num. obs.                  280560            280560        
## ===========================================================
## *** p < 0.001; ** p < 0.01; * p < 0.05 (or 0 outside the confidence interval).
```

```r
## SAVE SERIALIZED
fits.file <- sprintf('fit_compare_%s_pd%s_R%s_%s_d%s.rds', firm_i, nPeriods, R, m_x, d)
saveRDS(fit4mc, file=file.path(data_dir,fits.file))
## SAVE FORMATTED REGRESSION TABLE
html.file <- sprintf('%s_compare_tergm_results_pd%s_R%s_%s_d%s.html',  firm_i, nPeriods, R, m_x, d)
htmlreg(fit4mc, digits = 2, file=file.path(data_dir,html.file))
```

```
## The table was written to the file 'C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/awareness-cues-amj/qualtrics_compare_tergm_results_pd11_R100_m4_d2.html'.
```

```r
cat('finished successfully.')
```

```
## finished successfully.
```

