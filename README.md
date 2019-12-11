# Awareness Cues Computations

This repository contains R scripts and links to serialized data files for computations used in 

 * Downing, Kang, & Markman (Forthcoming) "What you don't see can hurt you: Awareness cues to profile indirect competitors." _Academy of Management Journal_. [https://doi.org/10.5465/amj.2018.0048](https://doi.org/10.5465/amj.2018.0048 "AMJ Article Link").


#### Acknowledgement
This research was supported in part by the Taiwan Ministry of Science and Technology: MOST 107-2410-H-009-052, MOST-105-2420-H-009-012-DR.

#### Background 
- [Intro to Network Analysis in R](https://statnet.org/trac/raw-attachment/wiki/Resources/introToSNAinR_sunbelt_2012_tutorial.pdf  "link")    
- [The R network Package](https://www.jstatsoft.org/index.php/jss/article/view/v024i02/v24i02.pdf  "download")    
- [Intro to Exponential Random Graph Models (ERGM)](http://ranger.uta.edu/~chqding/cse5301/classPapers/ExponentialRandomGraph.pdf  "link")    
- [Computing ERGMs in R](https://www.jstatsoft.org/index.php/jss/article/view/v024i03/v24i03.pdf  "download")
- [Bootstrapped ERGMs for Big Networks in R](https://arxiv.org/pdf/1708.02598.pdf  "link")
- [Temporal ERGMs (TERGM) in R](https://www.jstatsoft.org/index.php/jss/article/view/v083i06/v83i06.pdf "download")
- [Compare TERGM vs. SAOM](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/85F5AB6B84B2D49B05A5D24F4A018E88/S2050124218000267a.pdf/theoretical_and_empirical_comparison_of_the_temporal_exponential_random_graph_model_and_the_stochastic_actororiented_model.pdf "link")

#### Contents
- [Part 1: Prepare Data](#part-1-prepare-data  "Part 1")
- [Part 2: Create Global Competition Network](#part-2-create-global-competition-network  "Part 2")
- [Part 3: Compute Awareness Cues and Covariate Lists](#part-3-compute-awareness-cues-and-covariate-lists  "Part 3")
- [Part 4: Estimate TERGM](#part-4-estimate-tergm  "Part 4")
- [Part 5: Goodness of Fit](#part-5-goodness-of-fit  "Part 5")
- [Part 6: Estimation Algorithm](#part-6-estimation-algorithm  "Part 6")




# Part 1: Prepare Data

> [Back to Contents](#contents  "Back")


```r
##===============================
## SET YOUR DIRECTORIES:
##   This is the path to the folder where you saved the data file.
##   If you are using a Windows PC, use double backslash path separators "..\\dir\\subdir\\.."
##-------------------------------
## DIRECTORIES
data_dir <- "/set/data/dir"
work_dir <- "/set/working/dir"
img_dir  <- "/set/image/dir"
version_dir <- "/set/version/dir"
net_dir <- "/set/networks/dir"
sup_data_dir <- "/set/supplementalData/dir"
out_dir <- "/set/output/dir"

## set woring dir
setwd(work_dir)
```

Load `R` scripts:
- `amj_awareness_functions.R` loads functions for data processing; cached in environment as [list] `aaf`
- `amj_cb_data_prep.R` loads CrunchBase data tables; cached in environment as [list] `cb`
- `amj_sdc_coop.R` loads SDC data tables; cached in environment as [list] `sdc`
- `amj_firm_size_controls.R` loads firm size controls tables; cached in environment as [list] `si`
- `amj_institutional_holdings.R` loads institutional holdings tables; cached in environment as [list] `ih`
- `amj_make_full_graph.R` loads global competition network object; cached in environment as [list] `g.full`

This will take several minutes to complete.

```r
## LOAD DATA AND DEPENDENCIES
aaf    <- source(file.path(version_dir,'amj_awareness_functions.R'))$value    ## aaf: awareness functions
cb     <- source(file.path(version_dir,'amj_cb_data_prep.R'))$value           ## cb : CrunchBase
sdc    <- source(file.path(version_dir,'amj_sdc_coop.R'))$value               ## sdc: Thompson SDC
si     <- source(file.path(version_dir,'amj_firm_size_controls.R'))$value     ## si : size controls from mergent intellect
ih     <- source(file.path(version_dir,'amj_institutional_holdings.R'))$value ## ih : institutional holdings
g.full <- source(file.path(version_dir,'amj_make_full_graph.R'))$value        ## g.full : full competition graph
```

```
## 
## loading dataframes...done.
## cleaning data...
## Warning: All formats failed to parse. No formats found.
## Warning: All formats failed to parse. No formats found.
## Warning: All formats failed to parse. No formats found.
## 
  |                                                                       
  |===============                                                  |  24%
  |                                                                       
  |=================================================================| 100%
## reshaping acquisitions dataframe...                                                                  
  |                                                                       
  |=                                                                |   1% 
  |=================================================================| 100%
## done.
## clearing environment...
## 
## loading SDC cooperative relations data...
## loading firm size controls...
## loading institutional holdings data...
##
```

```r
print(summary(aaf))
```

```
##                       Length Class  Mode    
## makeGraph             1      -none- function
## generalistIndex       1      -none- function
## jobsToBeDone          1      -none- function
## mmcMarketsDf          1      -none- function
## mmcfromMarketConcat   1      -none- function
## coopConcatDf          1      -none- function
## coopFromConcat        1      -none- function
## .cov.coop             1      -none- function
## .cov.coopPast         1      -none- function
## .cov.age              1      -none- function
## .cov.mmc              1      -none- function
## .cov.dist             1      -none- function
## .cov.ipo              1      -none- function
## .cov.constraint       1      -none- function
## .cov.similarity       1      -none- function
## .cov.centrality       1      -none- function
## .cov.generalistIndex  1      -none- function
## setCovariates         1      -none- function
## nodeCollapseGraph     1      -none- function
## makePdNetwork         1      -none- function
## makePdGraph           1      -none- function
## getNetEcount          1      -none- function
## plotCompNetColPredict 1      -none- function
```

```r
print(summary(cb))
```

```
##                 Length Class      Mode    
## uuid             1     -none-     function
## match            1     -none-     function
## parseNonNa       1     -none-     function
## falsy            1     -none-     function
## relationBeganOn  1     -none-     function
## relationEndedOn  1     -none-     function
## readCsv          1     -none-     function
## fixDateYMD       1     -none-     function
## csv             13     -none-     list    
## co              33     data.frame list    
## co_comp         15     data.frame list    
## co_cust          7     data.frame list    
## co_parent        7     data.frame list    
## co_prod         12     data.frame list    
## co_acq          23     data.frame list    
## co_br           19     data.frame list    
## co_rou          21     data.frame list    
## co_ipo          20     data.frame list    
## fund             9     data.frame list    
## inv             14     data.frame list    
## inv_rou          3     data.frame list    
## inv_part         3     data.frame list
```





# Part 2: Create Global Competition Network

> [Back to Contents](#contents  "Back")


Create the full competition network (graph) for all competitive relations at all times. The following step after this will then create competition network panels with one competition network per time period by removing the relations and firms that didn't exist during that period.


```r
##==================================================
##
##  Make Full Graph
##
##--------------------------------------------------

cat('\nmaking full graph...')
```

```
## 
## making full graph...
```

```r
max.year <- 2016

## delete edges at or later than this date (the year after max.year)
exclude.date <- sprintf('%d-01-01', max.year+1)

## make graph
g.full <- aaf$makeGraph(comp = cb$co_comp, vertdf = cb$co)

## cut out confirmed dates >= 2016
g.full <- igraph::induced.subgraph(g.full, vids=V(g.full)[which(V(g.full)$founded_year <= max.year
                                                                | is.na(V(g.full)$founded_year)
                                                                | V(g.full)$founded_year=='' ) ] )
g.full <- igraph::delete.edges(g.full, E(g.full)[which(E(g.full)$relation_created_at >= exclude.date)])

## SIMPLIFY
g.full <- igraph::simplify(g.full, remove.loops=T,remove.multiple=T,
                           edge.attr.comb = list(weight='sum',
                                                 relation_began_on='max',
                                                 relation_ended_on='min'))

## save graph file
igraph::write.graph(graph = g.full, file=file.path(data_dir, "g_full.graphml"), format = 'graphml')

cat('done.\n')
```

```
## done.
```




# Part 3: Compute Awareness Cues and Covariate Lists

> [Back to Contents](#contents  "Back")

The main data preparation step involves using the full competition network and running a 3-step procedure to compute temporal panel data of competition networks and covariate arrays for each period.

The steps apply functions loaded in the `aaf` object for each period in the analysis time frame:
1. `aaf$nodeCollapseGraph(...)` Process acquisitions by transferring competitive relations from acquisition target to acquiring firm for each period
2. `aaf$makePdNetwork(...)` Filter the competitive relations and firms that existed within each period
3. `aaf$setCovariates(...)`  Compute node and edge covariates from the updated period competition network and set the covariates in this period's `network` object


```r
# ## set firms to create networks (focal firm or replication study focal firms)
firms.todo <- c(
  'qualtrics',
  ## other firm names to process
)

## -- settings --
d <- 2
yrpd <- 1
startYr <- 2005
endYr <- 2017            ## dropping first for memory term; actual dates 2007-2016
lg.cutoff <- 1100        ## large network size cutoff to save periods seprately 
force.overwrite <- FALSE ## if network files in directory should be overwritten
## --------------  


##
## run main network period creation loop
##
for (i in 1:length(firms.todo)) {
  if (i > 1) break ## only run first firm then stop

  name_i <- firms.todo[i]
  cat(sprintf('\n\n------------ %s -------------\n\n',name_i))
  periods <- seq(startYr,endYr,yrpd)
  company.name <- 'company_name_unique'
  g.base <- g.full  
  
  ## focal firm ego network sample
  g.d.sub <- igraph::make_ego_graph(graph = g.base, nodes = V(g.base)[V(g.base)$name==name_i], order = d, mode = 'all')[[1]]
  
  ## convert to network object
  net.d.sub <- asNetwork(g.d.sub)
  net <- net.d.sub
  net %n% 'ego' <- name_i
  
  ##------------------------------------------------------
  ##-------preprocess parent-subsidiary relationships-----
  ##----------node collapse like acquisitions-------------
  ##------------------------------------------------------
  cat(' collapsing parent-subsidiary relations...')
  ## load in manually checked parent-subsidiary relations
  # dfpar.xl <- file.path(sup_data_dir, 'private_firm_financials_all_firm_networks_delete_nodes.xlsx')
  # dfpar <- read_excel(dfpar.xl,  na = c('','-',"'-"))
  dfpar.csv <- file.path(sup_data_dir, 'private_firm_financials_all_firm_networks_delete_nodes.csv')
  dfpar <- read.csv(dfpar.csv, na.strings = c('','-',"'-","--"), stringsAsFactors = F)
  dfpar <- dfpar[!is.na(dfpar$parent), ]
  dfpar$parent <- str_to_lower(dfpar$parent)
  ## merge in parent uuid
  .parent.uuid <- cb$co[,c('company_name_unique','company_uuid')]
  names(.parent.uuid) <- c('parent_name_unique', 'parent_uuid') 
  dfpar <- merge(dfpar[,c('parent','firm')], .parent.uuid, by.x='parent', by.y='parent_name_unique', all.x=T, all.y=F)
  ## merge in firm uuid
  .firm.uuid <- cb$co[,c('company_name_unique','company_uuid')]
  dfpar <- merge(dfpar, .firm.uuid, by.x='firm', by.y='company_name_unique', all.x=T, all.y=F)
  ## CrunchBase parent relationships to collapse
  cb.par <- cb$co_parent[,c('company_name_unique','parent_name_unique','parent_org_uuid','org_uuid')]
  names(cb.par) <- c('firm','parent','parent_uuid','company_uuid')
  ## combine CrunchBase and manual parent-child mappings
  par.chi <- rbind(dfpar, cb.par)
  ## filter both parent, child in full graph
  gfuuid <- V(g.d.sub)$company_uuid
  par.chi <- par.chi[which(!is.na(par.chi$parent_uuid) & par.chi$parent_uuid %in% gfuuid 
                           & !is.na(par.chi$company_uuid) & par.chi$company_uuid %in% gfuuid), ]
  ## merge in founded_on date of parent
  par.chi <- merge(par.chi, cb$co[,c('company_name_unique','founded_on')], by.x='parent', by.y='company_name_unique', all.x=T, all.y=F)
  ## quasi-node collapse subsidiaries to parent nodes
  par.chi.nc <- data.frame(
    acquirer_uuid=par.chi$parent_uuid,
    acquiree_uuid=par.chi$company_uuid,
    acquired_on=par.chi$founded_on,   ## for parent-subsidiary mapping, just use parent company founded_on date
    stringsAsFactors = F
  )
  g.d.sub <- aaf$nodeCollapseGraph(g.d.sub, par.chi.nc, remove.isolates=T, verbose = T)
  
  cat('done.\n')
  
  ##_-----------------------------------------------------
  ##-------process pre-start-year acquisitions------------
  ##------------------------------------------------------
  acqs.pd <- cb$co_acq[cb$co_acq$acquired_on <= sprintf('%d-12-31',startYr-1), ]
  g.d.sub <- aaf$nodeCollapseGraph(g.d.sub, acqs.pd, remove.isolates=T, verbose = T)
  net.d.sub <- asNetwork(g.d.sub)
  cat(sprintf('v = %d, e = %d\n',vcount(g.d.sub),ecount(g.d.sub)))
  
  # ## subset to firms with employees count > 10
  # idx.employee <- which( !(V(g.d.sub)$employee_count %in% c('NA','-','1-10')) )
  # g.d.sub <- igraph::induced.subgraph(g.d.sub, vids = V(g.d.sub)[idx.employee])
  # cat(sprintf('filtered >10 employee count: v = %d, e = %d\n',vcount(g.d.sub),ecount(g.d.sub)))
  
  ##------------Network Time Period List--------------------
  nl <- list()
  
  for (t in 2:length(periods)) 
  {
    ## period dates
    cat(sprintf('\n\nmaking period %s-%s:\n', periods[t-1],periods[t]))
    t1 <- sprintf('%d-01-01',periods[t-1]) ## inclusive start date 'YYYY-MM-DD'
    t2 <- sprintf('%d-12-31',periods[t-1]) ## inclusive end date 'YYYY-MM-DD'
    
    ## check if period network file exists (skip if not force overwrite)
    file.rds <- file.path(net_dir,sprintf('%s_d%d_y%s.rds',name_i,d,periods[t-1]))
    if (!force.overwrite & file.exists(file.rds)) {
      cat(sprintf('file exists: %s\nskipping.\n', file.rds))
      next
    }
    
    ## period years indicate [start, end) -- start inclusive; end exclusive 
    
    ## 1. Node Collapse acquisitions within period
    acqs.pd <- cb$co_acq[cb$co_acq$acquired_on >= t1 & cb$co_acq$acquired_on <= t2, ]
    g.d.sub <- aaf$nodeCollapseGraph(g.d.sub, acqs.pd, verbose = T)
    
    ## 2. Subset Period Network
    nl[[t]] <- aaf$makePdNetwork(asNetwork(g.d.sub), periods[t-1], periods[t], isolates.remove = F) 
    
    ## 3. Set Covariates for updated Period Network
    nl[[t]] <- aaf$setCovariates(nl[[t]], periods[t-1], periods[t],
                                 acq=cb$co_acq, br=cb$co_br, ipo=cb$co_ipo, 
                                 rou=cb$co_rou, inv_rou=cb$inv_rou, inv=cb$inv,
                                 coop=sdc, ih=ih, size=si)
    
    ## save each period if large network (would exceed memory as full list of time periods)
    if (vcount(g.d.sub) >= lg.cutoff) {
      saveRDS(nl[[t]], file = file.rds)
      nv <- length(nl[[t]]$val)
      names(nv)[1] <- as.character(periods[t-1])
      write.csv(nv, file = file.path(out_dir, sprintf('%s_d%s.csv',name_i,d)),append = TRUE)
      nl[[t]] <- NULL ## remove from memory
    }
    
  }
  
  ##---------Small Networks: clean period list and save whole -----------
  if (vcount(g.d.sub) < lg.cutoff) 
  {
    ## ----drop null and skipped periods----
    nl.bak <- nl
    nl <- nl[which(sapply(nl, length)>0)]
    
    if (length(nl) > 1) {
      names(nl) <- periods[2:length(periods)]
    }

    ##--------------- GET TERGM NETS LIST -----------
    ## only nets with edges > 0
    if (length(nl) > 1) {
      nets.all <- nl[2:length(nl)]
    } else {
      nets.all <- nl
    }
    nets <- nets.all[ which(sapply(nets.all, aaf$getNetEcount) > 0) ]
    ## record network sizes
    write.csv(sapply(nets,function(x)length(x$val)), file = file.path(out_dir, sprintf('%s_d%s.csv',name_i,d)))
    
    #-------------------------------------------------
    
    ## CAREFUL TO OVERWRITE 
    saveRDS(nets, file = file.path(out_dir, sprintf('%s_d%d.rds',name_i,d)))
    
    ## plot covariate summary figures
    aaf$covSummaryPlot(nets, name_i, net_dir)
    
  }
  
}
```

```
## 
## 
## ------------ qualtrics -------------
## 
##  collapsing parent-subsidiary relations...processing acquisitions: 9 ...
## simplifying edges...done.
## simplifying edges...done.
## simplifying edges...done.
## simplifying edges...done.
## simplifying edges...done.
## simplifying edges...done.
## done.
## done.
## processing acquisitions: 0 ...
## done.
## v = 168, e = 380
## 
## 
## making period 2005-2006:
## processing acquisitions: 0 ...
## done.
## collecting edges and vertices to remove...done.
## computing age ...done
## computing multi-market contact (branch geographic overlap)...no firm branches in period. creating empty mmc matrix...done.
## done
## computing distances lag contact...done
## computing IPO status contact...done
## computing constraint...done
## computing inv.log.w.similarity...done
## computing centralities...
```

```
## done
## computing Generalist (vs Specialist) Index...
```

```
## Warning in igraph::edge.betweenness.community(g): At community.c:
## 460 :Membership vector will be selected based on the lowest modularity
## score.
```

```
## Warning in igraph::edge.betweenness.community(g): At community.c:
## 467 :Modularity calculation with weighted edge betweenness community
## detection might not make sense -- modularity treats edge weights as
## similarities while edge betwenness treats them as distances
```

```
## done
## computing Cooperative relations (alliance/JV)...concatenating current cooperative relations...
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |=================================================================| 100%
## done.
## computing current cooperative relations outerproduct matrix...done.
## concatenating past cooperative relations...
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |=================================================================| 100%
## done.
## computing past cooperative relations outerproduct matrix...done.
## done
## computing Employees...done
## computing Sales...done
## computing Category Cosine Similarity ...done
## computing Shared Competitor  ...done
## computing Shared Investors...
##  computing public firms common investors... done.
##  computing private firms common investors... done.
##  skipping mixed dyad public-private off-diagonal blocks.done
## 
## 
## making period 2006-2007:
## processing acquisitions: 0 ...
## done.
## collecting edges and vertices to remove...done.
## computing age ...done
## computing multi-market contact (branch geographic overlap)...no firm branches in period. creating empty mmc matrix...done.
## done
## computing distances lag contact...done
## computing IPO status contact...done
## computing constraint...done
## computing inv.log.w.similarity...done
## computing centralities...
```

```
## done
## computing Generalist (vs Specialist) Index...
```

```
## Warning in igraph::edge.betweenness.community(g): At community.c:
## 460 :Membership vector will be selected based on the lowest modularity
## score.
```

```
## Warning in igraph::edge.betweenness.community(g): At community.c:
## 467 :Modularity calculation with weighted edge betweenness community
## detection might not make sense -- modularity treats edge weights as
## similarities while edge betwenness treats them as distances
```

```
## done
## computing Cooperative relations (alliance/JV)...concatenating current cooperative relations...
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |=================================================================| 100%
## done.
## computing current cooperative relations outerproduct matrix...done.
## concatenating past cooperative relations...
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |=================================================================| 100%
## done.
## computing past cooperative relations outerproduct matrix...done.
## done
## computing Employees...done
## computing Sales...done
## computing Category Cosine Similarity ...done
## computing Shared Competitor  ...done
## computing Shared Investors...
##  computing public firms common investors... done.
##  computing private firms common investors... done.
##  skipping mixed dyad public-private off-diagonal blocks.done
## 
## 
## making period 2007-2008:
## processing acquisitions: 0 ...
## done.
## collecting edges and vertices to remove...done.
## computing age ...done
## computing multi-market contact (branch geographic overlap)...no firm branches in period. creating empty mmc matrix...done.
## done
## computing distances lag contact...done
## computing IPO status contact...done
## computing constraint...done
## computing inv.log.w.similarity...done
## computing centralities...
```

```
## done
## computing Generalist (vs Specialist) Index...
```

```
## Warning in igraph::edge.betweenness.community(g): At community.c:
## 460 :Membership vector will be selected based on the lowest modularity
## score.
```

```
## Warning in igraph::edge.betweenness.community(g): At community.c:
## 467 :Modularity calculation with weighted edge betweenness community
## detection might not make sense -- modularity treats edge weights as
## similarities while edge betwenness treats them as distances
```

```
## done
## computing Cooperative relations (alliance/JV)...concatenating current cooperative relations...
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |=================================================================| 100%
## done.
## computing current cooperative relations outerproduct matrix...done.
## concatenating past cooperative relations...
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |=================================================================| 100%
## done.
## computing past cooperative relations outerproduct matrix...done.
## done
## computing Employees...done
## computing Sales...done
## computing Category Cosine Similarity ...done
## computing Shared Competitor  ...done
## computing Shared Investors...
##  computing public firms common investors... done.
##  computing private firms common investors... done.
##  skipping mixed dyad public-private off-diagonal blocks.done
## 
## 
## making period 2008-2009:
## processing acquisitions: 1 ...
## simplifying edges...done.
## done.
## collecting edges and vertices to remove...done.
## computing age ...done
## computing multi-market contact (branch geographic overlap)...concatenating firm branch markets...
## 
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |======================                                           |  33%
  |                                                                       
  |===========================================                      |  67%
  |                                                                       
  |=================================================================| 100%
## computing MMC outer product matrix...done.
## done
## computing distances lag contact...done
## computing IPO status contact...done
## computing constraint...done
## computing inv.log.w.similarity...done
## computing centralities...
```

```
## done
## computing Generalist (vs Specialist) Index...
```

```
## Warning in igraph::edge.betweenness.community(g): At community.c:
## 460 :Membership vector will be selected based on the lowest modularity
## score.
```

```
## Warning in igraph::edge.betweenness.community(g): At community.c:
## 467 :Modularity calculation with weighted edge betweenness community
## detection might not make sense -- modularity treats edge weights as
## similarities while edge betwenness treats them as distances
```

```
## done
## computing Cooperative relations (alliance/JV)...concatenating current cooperative relations...
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |================================                                 |  50%
  |                                                                       
  |=================================================================| 100%
## done.
## computing current cooperative relations outerproduct matrix...done.
## concatenating past cooperative relations...
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |=================================================================| 100%
## done.
## computing past cooperative relations outerproduct matrix...done.
## done
## computing Employees...done
## computing Sales...done
## computing Category Cosine Similarity ...done
## computing Shared Competitor  ...done
## computing Shared Investors...
##  computing public firms common investors... done.
##  computing private firms common investors... done.
##  skipping mixed dyad public-private off-diagonal blocks.done
## 
## 
## making period 2009-2010:
## processing acquisitions: 0 ...
## done.
## collecting edges and vertices to remove...done.
## computing age ...done
## computing multi-market contact (branch geographic overlap)...concatenating firm branch markets...
## 
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |======================                                           |  33%
  |                                                                       
  |===========================================                      |  67%
  |                                                                       
  |=================================================================| 100%
## computing MMC outer product matrix...done.
## done
## computing distances lag contact...done
## computing IPO status contact...done
## computing constraint...done
## computing inv.log.w.similarity...done
## computing centralities...
```

```
## done
## computing Generalist (vs Specialist) Index...
```

```
## Warning in igraph::edge.betweenness.community(g): At community.c:
## 460 :Membership vector will be selected based on the lowest modularity
## score.
```

```
## Warning in igraph::edge.betweenness.community(g): At community.c:
## 467 :Modularity calculation with weighted edge betweenness community
## detection might not make sense -- modularity treats edge weights as
## similarities while edge betwenness treats them as distances
```

```
## done
## computing Cooperative relations (alliance/JV)...concatenating current cooperative relations...
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |================================                                 |  50%
  |                                                                       
  |=================================================================| 100%
## done.
## computing current cooperative relations outerproduct matrix...done.
## concatenating past cooperative relations...
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |================================                                 |  50%
  |                                                                       
  |=================================================================| 100%
## done.
## computing past cooperative relations outerproduct matrix...done.
## done
## computing Employees...done
## computing Sales...done
## computing Category Cosine Similarity ...done
## computing Shared Competitor  ...done
## computing Shared Investors...
##  computing public firms common investors... done.
##  computing private firms common investors... done.
##  skipping mixed dyad public-private off-diagonal blocks.done
## 
## 
## making period 2010-2011:
## processing acquisitions: 0 ...
## done.
## collecting edges and vertices to remove...done.
## computing age ...done
## computing multi-market contact (branch geographic overlap)...concatenating firm branch markets...
## 
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |======================                                           |  33%
  |                                                                       
  |===========================================                      |  67%
  |                                                                       
  |=================================================================| 100%
## computing MMC outer product matrix...done.
## done
## computing distances lag contact...done
## computing IPO status contact...done
## computing constraint...done
## computing inv.log.w.similarity...done
## computing centralities...
```

```
## done
## computing Generalist (vs Specialist) Index...
```

```
## Warning in igraph::edge.betweenness.community(g): At community.c:
## 460 :Membership vector will be selected based on the lowest modularity
## score.
```

```
## Warning in igraph::edge.betweenness.community(g): At community.c:
## 467 :Modularity calculation with weighted edge betweenness community
## detection might not make sense -- modularity treats edge weights as
## similarities while edge betwenness treats them as distances
```

```
## done
## computing Cooperative relations (alliance/JV)...concatenating current cooperative relations...
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |======================                                           |  33%
  |                                                                       
  |===========================================                      |  67%
  |                                                                       
  |=================================================================| 100%
## done.
## computing current cooperative relations outerproduct matrix...done.
## concatenating past cooperative relations...
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |================================                                 |  50%
  |                                                                       
  |=================================================================| 100%
## done.
## computing past cooperative relations outerproduct matrix...done.
## done
## computing Employees...done
## computing Sales...done
## computing Category Cosine Similarity ...done
## computing Shared Competitor  ...done
## computing Shared Investors...
##  computing public firms common investors... done.
##  computing private firms common investors... done.
##  skipping mixed dyad public-private off-diagonal blocks.done
## 
## 
## making period 2011-2012:
## processing acquisitions: 1 ...
## simplifying edges...done.
## done.
## collecting edges and vertices to remove...done.
## computing age ...done
## computing multi-market contact (branch geographic overlap)...concatenating firm branch markets...
## 
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |======================                                           |  33%
  |                                                                       
  |===========================================                      |  67%
  |                                                                       
  |=================================================================| 100%
## computing MMC outer product matrix...done.
## done
## computing distances lag contact...done
## computing IPO status contact...done
## computing constraint...done
## computing inv.log.w.similarity...done
## computing centralities...
```

```
## done
## computing Generalist (vs Specialist) Index...
```

```
## Warning in igraph::edge.betweenness.community(g): At community.c:
## 460 :Membership vector will be selected based on the lowest modularity
## score.
```

```
## Warning in igraph::edge.betweenness.community(g): At community.c:
## 467 :Modularity calculation with weighted edge betweenness community
## detection might not make sense -- modularity treats edge weights as
## similarities while edge betwenness treats them as distances
```

```
## done
## computing Cooperative relations (alliance/JV)...concatenating current cooperative relations...
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |======================                                           |  33%
  |                                                                       
  |===========================================                      |  67%
  |                                                                       
  |=================================================================| 100%
## done.
## computing current cooperative relations outerproduct matrix...done.
## concatenating past cooperative relations...
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |======================                                           |  33%
  |                                                                       
  |===========================================                      |  67%
  |                                                                       
  |=================================================================| 100%
## done.
## computing past cooperative relations outerproduct matrix...done.
## done
## computing Employees...done
## computing Sales...done
## computing Category Cosine Similarity ...done
## computing Shared Competitor  ...done
## computing Shared Investors...
##  computing public firms common investors... done.
##  computing private firms common investors... done.
##  skipping mixed dyad public-private off-diagonal blocks.done
## 
## 
## making period 2012-2013:
## processing acquisitions: 0 ...
## done.
## collecting edges and vertices to remove...done.
## computing age ...done
## computing multi-market contact (branch geographic overlap)...concatenating firm branch markets...
## 
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |======================                                           |  33%
  |                                                                       
  |===========================================                      |  67%
  |                                                                       
  |=================================================================| 100%
## computing MMC outer product matrix...done.
## done
## computing distances lag contact...done
## computing IPO status contact...done
## computing constraint...done
## computing inv.log.w.similarity...done
## computing centralities...
```

```
## done
## computing Generalist (vs Specialist) Index...
```

```
## Warning in igraph::edge.betweenness.community(g): At community.c:
## 460 :Membership vector will be selected based on the lowest modularity
## score.
```

```
## Warning in igraph::edge.betweenness.community(g): At community.c:
## 467 :Modularity calculation with weighted edge betweenness community
## detection might not make sense -- modularity treats edge weights as
## similarities while edge betwenness treats them as distances
```

```
## done
## computing Cooperative relations (alliance/JV)...concatenating current cooperative relations...
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |======================                                           |  33%
  |                                                                       
  |===========================================                      |  67%
  |                                                                       
  |=================================================================| 100%
## done.
## computing current cooperative relations outerproduct matrix...done.
## concatenating past cooperative relations...
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |======================                                           |  33%
  |                                                                       
  |===========================================                      |  67%
  |                                                                       
  |=================================================================| 100%
## done.
## computing past cooperative relations outerproduct matrix...done.
## done
## computing Employees...done
## computing Sales...done
## computing Category Cosine Similarity ...done
## computing Shared Competitor  ...done
## computing Shared Investors...
##  computing public firms common investors... done.
##  computing private firms common investors... done.
##  skipping mixed dyad public-private off-diagonal blocks.done
## 
## 
## making period 2013-2014:
## processing acquisitions: 0 ...
## done.
## collecting edges and vertices to remove...done.
## computing age ...done
## computing multi-market contact (branch geographic overlap)...concatenating firm branch markets...
## 
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |======================                                           |  33%
  |                                                                       
  |===========================================                      |  67%
  |                                                                       
  |=================================================================| 100%
## computing MMC outer product matrix...done.
## done
## computing distances lag contact...done
## computing IPO status contact...done
## computing constraint...done
## computing inv.log.w.similarity...done
## computing centralities...
```

```
## done
## computing Generalist (vs Specialist) Index...
```

```
## Warning in igraph::edge.betweenness.community(g): At community.c:
## 460 :Membership vector will be selected based on the lowest modularity
## score.
```

```
## Warning in igraph::edge.betweenness.community(g): At community.c:
## 467 :Modularity calculation with weighted edge betweenness community
## detection might not make sense -- modularity treats edge weights as
## similarities while edge betwenness treats them as distances
```

```
## done
## computing Cooperative relations (alliance/JV)...concatenating current cooperative relations...
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |======================                                           |  33%
  |                                                                       
  |===========================================                      |  67%
  |                                                                       
  |=================================================================| 100%
## done.
## computing current cooperative relations outerproduct matrix...done.
## concatenating past cooperative relations...
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |======================                                           |  33%
  |                                                                       
  |===========================================                      |  67%
  |                                                                       
  |=================================================================| 100%
## done.
## computing past cooperative relations outerproduct matrix...done.
## done
## computing Employees...done
## computing Sales...done
## computing Category Cosine Similarity ...done
## computing Shared Competitor  ...done
## computing Shared Investors...
##  computing public firms common investors... done.
##  computing private firms common investors... done.
##  skipping mixed dyad public-private off-diagonal blocks.done
## 
## 
## making period 2014-2015:
## processing acquisitions: 0 ...
## done.
## collecting edges and vertices to remove...done.
## computing age ...done
## computing multi-market contact (branch geographic overlap)...concatenating firm branch markets...
## 
  |                                                                       
  |==                                                               |   3%
  |                                                                       
  |=================================================================| 100%
## computing MMC outer product matrix...done.
## done
## computing distances lag contact...done
## computing IPO status contact...done
## computing constraint...done
## computing inv.log.w.similarity...done
## computing centralities...
```

```
## done
## computing Generalist (vs Specialist) Index...
```

```
## Warning in igraph::edge.betweenness.community(g): At community.c:
## 460 :Membership vector will be selected based on the lowest modularity
## score.
```

```
## Warning in igraph::edge.betweenness.community(g): At community.c:
## 467 :Modularity calculation with weighted edge betweenness community
## detection might not make sense -- modularity treats edge weights as
## similarities while edge betwenness treats them as distances
```

```
## done
## computing Cooperative relations (alliance/JV)...concatenating current cooperative relations...
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |======================                                           |  33%
  |                                                                       
  |===========================================                      |  67%
  |                                                                       
  |=================================================================| 100%
## done.
## computing current cooperative relations outerproduct matrix...done.
## concatenating past cooperative relations...
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |======================                                           |  33%
  |                                                                       
  |===========================================                      |  67%
  |                                                                       
  |=================================================================| 100%
## done.
## computing past cooperative relations outerproduct matrix...done.
## done
## computing Employees...done
## computing Sales...done
## computing Category Cosine Similarity ...done
## computing Shared Competitor  ...done
## computing Shared Investors...
##  computing public firms common investors... done.
##  computing private firms common investors... done.
##  skipping mixed dyad public-private off-diagonal blocks.done
## 
## 
## making period 2015-2016:
## processing acquisitions: 0 ...
## done.
## collecting edges and vertices to remove...done.
## computing age ...done
## computing multi-market contact (branch geographic overlap)...concatenating firm branch markets...
## 
  |                                                                       
  |===                                                              |   5%
  |                                                                       
  |=================================================================| 100%
## computing MMC outer product matrix...done.
## done
## computing distances lag contact...done
## computing IPO status contact...done
## computing constraint...done
## computing inv.log.w.similarity...done
## computing centralities...
```

```
## done
## computing Generalist (vs Specialist) Index...
```

```
## Warning in igraph::edge.betweenness.community(g): At community.c:
## 460 :Membership vector will be selected based on the lowest modularity
## score.
```

```
## Warning in igraph::edge.betweenness.community(g): At community.c:
## 467 :Modularity calculation with weighted edge betweenness community
## detection might not make sense -- modularity treats edge weights as
## similarities while edge betwenness treats them as distances
```

```
## done
## computing Cooperative relations (alliance/JV)...concatenating current cooperative relations...
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |======================                                           |  33%
  |                                                                       
  |===========================================                      |  67%
  |                                                                       
  |=================================================================| 100%
## done.
## computing current cooperative relations outerproduct matrix...done.
## concatenating past cooperative relations...
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |======================                                           |  33%
  |                                                                       
  |===========================================                      |  67%
  |                                                                       
  |=================================================================| 100%
## done.
## computing past cooperative relations outerproduct matrix...done.
## done
## computing Employees...done
## computing Sales...done
## computing Category Cosine Similarity ...done
## computing Shared Competitor  ...done
## computing Shared Investors...
##  computing public firms common investors... done.
##  computing private firms common investors... done.
##  skipping mixed dyad public-private off-diagonal blocks.done
## 
## 
## making period 2016-2017:
## processing acquisitions: 3 ...
## simplifying edges...done.
## simplifying edges...done.
## simplifying edges...done.
## done.
## collecting edges and vertices to remove...done.
## computing age ...done
## computing multi-market contact (branch geographic overlap)...concatenating firm branch markets...
## 
  |                                                                       
  |==                                                               |   3%
  |                                                                       
  |=================================================================| 100%
## computing MMC outer product matrix...done.
## done
## computing distances lag contact...done
## computing IPO status contact...done
## computing constraint...done
## computing inv.log.w.similarity...done
## computing centralities...
```

```
## done
## computing Generalist (vs Specialist) Index...
```

```
## Warning in igraph::edge.betweenness.community(g): At community.c:
## 460 :Membership vector will be selected based on the lowest modularity
## score.
```

```
## Warning in igraph::edge.betweenness.community(g): At community.c:
## 467 :Modularity calculation with weighted edge betweenness community
## detection might not make sense -- modularity treats edge weights as
## similarities while edge betwenness treats them as distances
```

```
## done
## computing Cooperative relations (alliance/JV)...concatenating current cooperative relations...
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |======================                                           |  33%
  |                                                                       
  |===========================================                      |  67%
  |                                                                       
  |=================================================================| 100%
## done.
## computing current cooperative relations outerproduct matrix...done.
## concatenating past cooperative relations...
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |======================                                           |  33%
  |                                                                       
  |===========================================                      |  67%
  |                                                                       
  |=================================================================| 100%
## done.
## computing past cooperative relations outerproduct matrix...done.
## done
## computing Employees...done
## computing Sales...done
## computing Category Cosine Similarity ...done
## computing Shared Competitor  ...done
## computing Shared Investors...
##  computing public firms common investors... done.
##  computing private firms common investors... done.
##  skipping mixed dyad public-private off-diagonal blocks.done
```

```r
## ----drop null and skipped periods----
nl.bak <- nl
nl <- nl[which(sapply(nl, length)>0)]

if (length(nl) > 1) {
  names(nl) <- periods[2:length(periods)]
}

##--------------- GET TERGM NETS LIST -----------
## only nets with edges > 0
if (length(nl) > 1) {
  nets.all <- nl[2:length(nl)]
} else {
  nets.all <- nl
}
nets <- nets.all[ which(sapply(nets.all, aaf$getNetEcount) > 0) ]
## record network sizes
write.csv(sapply(nets,function(x)length(x$val)), file = file.path(data_dir,sprintf('%s_d%s.csv',name_i,d)))

#-------------------------------------------------

## Save serialized data file of all networks and covariates lists 
file.rds <- file.path(data_dir,sprintf('%s_d%d.rds',name_i,d))
saveRDS(nets, file = file.rds)
```

Next, in [Part 4](#part-4-estimate-tergm  "Part 4") compute a TERGM with the data lists that were just computed.


> [Back to Contents](#contents  "Back")








# Part 4: Estimate TERGM

> [Back to Contents](#contents  "Back")

In this repository's `R` directory, refer to the R script `amj_TERGM_estimate.R`. 

Download and save the following RDS (serialized) data file     
- [Firm Network Serialized Data](https://drive.google.com/drive/folders/18bE99uD-TVgkeGR501Sh_VvAhc6D6Bgn?usp=sharing "Firm Network Serialized Data")

in the same directory that you save the above script. You can run the script in its entirety to get the results.

Set the name of the directory where you saved the data file:

```r
##===============================
## SET YOUR DATA DIRECTORY:
##   This is the path to the folder where you saved the data file.
##   If you are using a Windows PC, use double backslash path separators "..\\dir\\subdir\\.."
##-------------------------------
data_dir <- '/set/your/data/directory/here'
```

Set parameters for the analysis and load the data into memory:

```r
## analysis parameters
firm_i <- 'qualtrics'  ## focal firm
d <- 2                 ## ego network theshold (order)

## load RDS data file into memory as a list of networks
# data_file <- file.path(data_dir,sprintf('%s_d%s.rds',firm_i,d))
data_file <- "data_file_name.rds"  ## set data file name
nets.all <- readRDS(data_file)
len <- length(nets.all)

## set number of time periods
nPeriods <- 8  ## any number from 3 to 11 is OK, but use 8 to compare results with example

## subset network periods
nets <- nets.all[(len-nPeriods+1):len]
```

Set the model formulas for the Part 1 tutorial:

```r
m0 <-   nets ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed=T) + 
  nodematch("ipo_status", diff = F) + 
  nodematch("state_code", diff = F) + 
  nodecov("age") + absdiff("age") + 
  memory(type = "stability", lag = 1)

m1 <-   nets ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed=T) + 
  nodematch("ipo_status", diff = F) + 
  nodematch("state_code", diff = F) + 
  nodecov("age") + absdiff("age") + 
  nodecov("cent_deg") +
  memory(type = "stability", lag = 1) + 
  nodecov("genidx_multilevel") + 
  nodecov("cent_pow_n0_4") + absdiff("cent_pow_n0_4") + 
  cycle(3) + cycle(4) 
```

Set the number of bootstrap replications. According to [Leifeld, Cranmer, & Desmarais (2018)](https://www.jstatsoft.org/article/view/v083i06 "Temporal Exponential Random Graph Models with btergm") :
- Roughly 100 is enough for an approximate estimate
- On the order of 1000 or more for reporting results


```r
R <- 100  ## enough for a rough estimate
```

Compute the first model `m0` and save to disk as an RDS (serialized) file:

This may take a while to run (a few minutes to a couple hours)  depending upon:
- network size (number of nodes per network)
- number of periods (number of network panels)
- complexity of change statistics to compute for the predictors in the model 

```r
## set pseudorandom number generator seed for reproducibility
set.seed(1111)
## estimate the TERGM with bootstrapped PMLE
fit0 <- btergm(get('m0'), R=R, parallel = "multicore", ncpus = detectCores())
```

```
## 
## Initial dimensions of the network and covariates:
##              t=2 t=3 t=4 t=5 t=6 t=7 t=8
## nets (row)   180 180 180 180 180 180 180
## nets (col)   180 180 180 180 180 180 180
## memory (row) 180 180 180 180 180 180 180
## memory (col) 180 180 180 180 180 180 180
## 
## All networks are conformable.
## 
## Dimensions of the network and covariates after adjustment:
##              t=2 t=3 t=4 t=5 t=6 t=7 t=8
## nets (row)   180 180 180 180 180 180 180
## nets (col)   180 180 180 180 180 180 180
## memory (row) 180 180 180 180 180 180 180
## memory (col) 180 180 180 180 180 180 180
## 
## Starting pseudolikelihood estimation with 100 bootstrapping replications using multicore forking on 4 cores...
## Done.
```

```r
## SAVE SERIALIZED DATA
fit0_file <- file.path(data_dir,sprintf('fit_%s_pd%s_R%s_%s.rds', firm_i, nPeriods, R, 'm0'))
saveRDS(fit0, file=fit0_file)
```

Compute the second model `m1` and save to disk as an RDS (serialized) file:

```r
## set pseudorandom number generator seed for reproducibility
set.seed(1111)
## estimate the TERGM with bootstrapped PMLE
fit1 <- btergm(get('m1'), R=R, parallel = "multicore", ncpus = detectCores())  
```

```
## 
## Initial dimensions of the network and covariates:
##              t=2 t=3 t=4 t=5 t=6 t=7 t=8
## nets (row)   180 180 180 180 180 180 180
## nets (col)   180 180 180 180 180 180 180
## memory (row) 180 180 180 180 180 180 180
## memory (col) 180 180 180 180 180 180 180
## 
## All networks are conformable.
## 
## Dimensions of the network and covariates after adjustment:
##              t=2 t=3 t=4 t=5 t=6 t=7 t=8
## nets (row)   180 180 180 180 180 180 180
## nets (col)   180 180 180 180 180 180 180
## memory (row) 180 180 180 180 180 180 180
## memory (col) 180 180 180 180 180 180 180
## 
## Starting pseudolikelihood estimation with 100 bootstrapping replications using multicore forking on 4 cores...
## Done.
```

```r
## SAVE SERIALIZED DATA
fit1_file <- file.path(data_dir,sprintf('fit_%s_pd%s_R%s_%s.rds', firm_i, nPeriods, R, 'm1'))
saveRDS(fit1, file=fit1_file)
```

Create a list of model fits. Print the regression table to screen and save it as a formatted HTML file.
You should see results like these:


```r
## Cache model fits list
fits <- list(Model_0=fit0,Model_1=fit1)

## Echo model comparison table to screen
screenreg(fits, digits = 3)
```

```
## 
## =============================================================
##                            Model_0           Model_1         
## -------------------------------------------------------------
## edges                         -0.180            -1.916 *     
##                            [-1.749;  1.067]  [-3.590; -0.606]
## gwesp.fixed.0                  0.813 *           0.449 *     
##                            [ 0.658;  0.908]  [ 0.152;  0.625]
## gwdegree                      -1.145 *          -0.504       
##                            [-1.910; -0.527]  [-1.454;  0.337]
## nodematch.ipo_status           0.132            -0.176       
##                            [-0.259;  1.070]  [-0.659;  0.777]
## nodematch.state_code          -0.593 *          -0.413 *     
##                            [-0.770; -0.387]  [-0.710; -0.048]
## nodecov.age                   -0.139 *          -0.119 *     
##                            [-0.199; -0.065]  [-0.159; -0.051]
## absdiff.age                    0.157 *           0.142 *     
##                            [ 0.085;  0.216]  [ 0.075;  0.177]
## edgecov.memory[[i]]            5.229 *           4.907 *     
##                            [ 4.908;  5.769]  [ 4.605;  5.482]
## nodecov.cent_deg                                 0.009       
##                                              [-0.071;  0.067]
## nodecov.genidx_multilevel                        1.526 *     
##                                              [ 0.401;  2.503]
## nodecov.cent_pow_n0_4                            0.008       
##                                              [-0.117;  0.176]
## absdiff.cent_pow_n0_4                            0.203 *     
##                                              [ 0.098;  0.450]
## cycle3                                           0.349 *     
##                                              [ 0.121;  0.615]
## cycle4                                          -0.010       
##                                              [-0.054;  0.035]
## -------------------------------------------------------------
## Num. obs.                  13787             67288           
## =============================================================
## * 0 outside the confidence interval
```

```r
## SAVE FORMATTED REGRESSION TABLE
compare_file <- file.path(data_dir,sprintf('%s_tergm_results_pd%s_R%s_%s.html', firm_i, nPeriods, R, 'm0-m1'))
htmlreg(fits, digits = 3, file=compare_file)

## 
## finished successfully.
```

# Part 5: Goodness of Fit

> [Back to Contents](#contents  "Back")

Continue using the R script `amj_TERGM_estimate.R`. 

Check goodness of fit for the following diagnostic statistics by simulating `nsim` number of random networks from model `m0`:
- `dsp` dyad-wise shared partners
- `esp` edge-wise shared partners
- `degree` degree distribution
- `geodesic` shortest path distribution

On the order of 1,000 to 10,000 networks should be sampled for reporting goodness of fit results, but for instruction we briefly simulate 30 per period.

```r
gof0 <- gof(fit0, statistics=c(dsp, esp, deg, geodesic), nsim=30)
```

```
## 
## Starting GOF assessment on a single computing core....
## 
## Initial dimensions of the network and covariates:

##              t=2 t=3 t=4 t=5 t=6 t=7 t=8
## nets (row)   180 180 180 180 180 180 180
## nets (col)   180 180 180 180 180 180 180
## memory (row) 180 180 180 180 180 180 180
## memory (col) 180 180 180 180 180 180 180

## 
## All networks are conformable.
## 
## Dimensions of the network and covariates after adjustment:

##              t=2 t=3 t=4 t=5 t=6 t=7 t=8
## nets (row)   180 180 180 180 180 180 180
## nets (col)   180 180 180 180 180 180 180
## memory (row) 180 180 180 180 180 180 180
## memory (col) 180 180 180 180 180 180 180

## 
## No 'target' network(s) provided. Using networks on the left-hand side of the model formula as observed networks.
## Simulating 30 networks from the following formula:
##  nets[[1]] ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + nodematch("ipo_status", diff = F) + nodematch("state_code", diff = F) + nodecov("age") + absdiff("age") + edgecov(memory[[1]])
## Simulating 30 networks from the following formula:
##  nets[[2]] ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + nodematch("ipo_status", diff = F) + nodematch("state_code", diff = F) + nodecov("age") + absdiff("age") + edgecov(memory[[2]])
## Simulating 30 networks from the following formula:
##  nets[[3]] ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + nodematch("ipo_status", diff = F) + nodematch("state_code", diff = F) + nodecov("age") + absdiff("age") + edgecov(memory[[3]])
## Simulating 30 networks from the following formula:
##  nets[[4]] ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + nodematch("ipo_status", diff = F) + nodematch("state_code", diff = F) + nodecov("age") + absdiff("age") + edgecov(memory[[4]])
## Simulating 30 networks from the following formula:
##  nets[[5]] ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + nodematch("ipo_status", diff = F) + nodematch("state_code", diff = F) + nodecov("age") + absdiff("age") + edgecov(memory[[5]])
## Simulating 30 networks from the following formula:
##  nets[[6]] ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + nodematch("ipo_status", diff = F) + nodematch("state_code", diff = F) + nodecov("age") + absdiff("age") + edgecov(memory[[6]])
## Simulating 30 networks from the following formula:
##  nets[[7]] ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + nodematch("ipo_status", diff = F) + nodematch("state_code", diff = F) + nodecov("age") + absdiff("age") + edgecov(memory[[7]])
## 7 networks from which simulations are drawn were provided.
## Processing statistic: Dyad-wise shared partners
## Processing statistic: Edge-wise shared partners
## Processing statistic: Degree
## Processing statistic: Geodesic distances
```

```r
print(gof0)
```

```
## Dyad-wise shared partners

##     obs: mean median   min   max  sim: mean median   min   max    Pr(>z)
## 0  2.6815e+04  26890 23734 30776 2.7695e+04  27529 23360 31484 0.4762126
## 1  4.1563e+03   4294   972  6534 3.4573e+03   3760   556  6896 0.4697068
## 2  9.1171e+02    698   366  1456 7.8336e+02    661   130  1766 0.5246196
## 3  1.9686e+02    198    72   298 1.7090e+02    141    28   374 0.5320781
## 4  8.2857e+01     90    26   130 6.5771e+01     73     4   156 0.3239972
## 5  3.1143e+01     28     4    58 2.5724e+01     26     0    78 0.5430098
## 6  9.1429e+00      8     0    16 8.3524e+00      8     0    28 0.7849971
## 7  6.2857e+00      6     0    14 4.2381e+00      2     0    16 0.3638171
## 8  2.5714e+00      2     0     6 1.8952e+00      0     0    10 0.5446789
## 9  5.7143e-01      0     0     4 1.0571e+00      0     0     6 0.4332566
## 10 3.1429e+00      2     0     8 2.2381e+00      0     0    10 0.5161993
## 11 8.5714e-01      0     0     2 8.1905e-01      0     0     6 0.9291793
## 12 2.8571e-01      0     0     2 2.5714e-01      0     0     4 0.9246306
## 13 0.0000e+00      0     0     0 1.9048e-02      0     0     2 0.1577941
## 14 8.5714e-01      0     0     2 4.1905e-01      0     0     2 0.3226778
## 15 0.0000e+00      0     0     0 1.5238e-01      0     0     2 4.804e-05
## 16 0.0000e+00      0     0     0 0.0000e+00      0     0     0 1.0000000
## 17 0.0000e+00      0     0     0 2.3810e-01      0     0     2 2.735e-07
## 18 2.8571e-01      0     0     2 4.7619e-02      0     0     2 0.4374178
## 19 2.8571e-01      0     0     2 2.5714e-01      0     0     2 0.9244243
## 20 0.0000e+00      0     0     0 3.1429e-01      0     0     2 2.366e-09
## 21 2.8571e-01      0     0     2 1.7143e-01      0     0     2 0.7050653
## 22 5.7143e-01      0     0     2 8.5714e-02      0     0     2 0.2366290
## 23 2.8571e-01      0     0     2 2.7619e-01      0     0     2 0.9747830
## 24 2.8571e-01      0     0     2 4.6667e-01      0     0     2 0.5560331
## 25 8.5714e-01      0     0     2 2.8571e-01      0     0     4 0.2088791
## 26 0.0000e+00      0     0     0 1.1429e-01      0     0     2 0.0004607
## 27 0.0000e+00      0     0     0 2.8571e-02      0     0     2 0.0832610
## 28 0.0000e+00      0     0     0 0.0000e+00      0     0     0 1.0000000
##       
## 0     
## 1     
## 2     
## 3     
## 4     
## 5     
## 6     
## 7     
## 8     
## 9     
## 10    
## 11    
## 12    
## 13    
## 14    
## 15 ***
## 16    
## 17 ***
## 18    
## 19    
## 20 ***
## 21    
## 22    
## 23    
## 24    
## 25    
## 26 ***
## 27 .  
## 28

## 
## Note: Small p-values indicate a significant difference 
##       between simulations and observed network(s).

## Edge-wise shared partners

##    obs: mean median min max  sim: mean median min max    Pr(>z)    
## 0  150.85714    166  92 194 150.190476    157  72 206 0.9648816    
## 1  167.14286    156 102 234 175.876190    178  64 300 0.7020008    
## 2   84.57143     78  42 126  78.219048     72  18 146 0.6350366    
## 3   62.85714     74  30  92  55.647619     49  18 120 0.4986991    
## 4   32.85714     36  14  54  28.095238     30   2  70 0.4291294    
## 5   18.57143     18   2  34  14.647619     13   0  48 0.4844099    
## 6    7.14286      8   0  12   6.047619      6   0  24 0.6165301    
## 7    4.57143      4   0  10   3.228571      2   0  12 0.4255438    
## 8    1.71429      2   0   4   1.247619      0   0   6 0.5221666    
## 9    0.57143      0   0   4   1.019048      0   0   6 0.4683639    
## 10   3.14286      2   0   8   2.238095      0   0  10 0.5161993    
## 11   0.57143      0   0   2   0.438095      0   0   6 0.7341107    
## 12   0.00000      0   0   0   0.085714      0   0   4 0.0063706 ** 
## 13   0.00000      0   0   0   0.000000      0   0   0 1.0000000    
## 14   0.00000      0   0   0   0.000000      0   0   0 1.0000000    
## 15   0.00000      0   0   0   0.019048      0   0   2 0.1577941    
## 16   0.00000      0   0   0   0.000000      0   0   0 1.0000000    
## 17   0.00000      0   0   0   0.238095      0   0   2 2.735e-07 ***
## 18   0.28571      0   0   2   0.047619      0   0   2 0.4374178    
## 19   0.28571      0   0   2   0.257143      0   0   2 0.9244243    
## 20   0.00000      0   0   0   0.314286      0   0   2 2.366e-09 ***
## 21   0.28571      0   0   2   0.171429      0   0   2 0.7050653    
## 22   0.57143      0   0   2   0.085714      0   0   2 0.2366290    
## 23   0.28571      0   0   2   0.276190      0   0   2 0.9747830    
## 24   0.28571      0   0   2   0.466667      0   0   2 0.5560331    
## 25   0.85714      0   0   2   0.285714      0   0   4 0.2088791    
## 26   0.00000      0   0   0   0.114286      0   0   2 0.0004607 ***
## 27   0.00000      0   0   0   0.028571      0   0   2 0.0832610 .  
## 28   0.00000      0   0   0   0.000000      0   0   0 1.0000000

## 
## Note: Small p-values indicate a significant difference 
##       between simulations and observed network(s).

## Degree

##    obs: mean median min max  sim: mean median min max    Pr(>z)    
## 0   50.71429     46  12 105 55.8238095     53  11 110 0.7432710    
## 1   45.00000     45  22  61 32.3761905     32  19  48 0.0810802 .  
## 2   21.28571     24  14  26 26.9380952     28  10  44 0.0316864 *  
## 3   18.28571     16  11  25 19.0714286     19   5  36 0.7410391    
## 4    9.42857      9   7  12 11.1190476     11   1  22 0.0706291 .  
## 5   11.42857     12   6  17 10.3476190     10   1  20 0.5105855    
## 6    8.57143     10   4  12  7.8333333      7   0  17 0.5384777    
## 7    2.00000      2   1   3  4.1047619      4   0  13 0.0001864 ***
## 8    3.71429      4   3   5  3.0380952      3   0   7 0.0572705 .  
## 9    1.42857      2   0   2  1.9142857      2   0   6 0.1594537    
## 10   1.00000      1   0   2  1.4095238      1   0   5 0.2407204    
## 11   1.42857      1   1   2  1.1380952      1   0   6 0.2142620    
## 12   0.42857      0   0   1  0.4190476      0   0   2 0.9644844    
## 13   0.42857      0   0   1  0.3000000      0   0   2 0.5526279    
## 14   0.28571      0   0   1  0.3285714      0   0   2 0.8268817    
## 15   0.28571      0   0   1  0.3952381      0   0   2 0.5798620    
## 16   0.00000      0   0   0  0.1523810      0   0   1 4.326e-09 ***
## 17   0.00000      0   0   0  0.1142857      0   0   2 1.426e-06 ***
## 18   0.28571      0   0   1  0.1095238      0   0   1 0.3783896    
## 19   0.28571      0   0   1  0.0666667      0   0   1 0.2809840    
## 20   0.14286      0   0   1  0.1285714      0   0   1 0.9244243    
## 21   0.14286      0   0   1  0.1428571      0   0   2 1.0000000    
## 22   0.57143      1   0   1  0.2904762      0   0   2 0.2166996    
## 23   0.00000      0   0   0  0.1333333      0   0   1 4.698e-08 ***
## 24   0.00000      0   0   0  0.0714286      0   0   1 8.464e-05 ***
## 25   0.14286      0   0   1  0.0714286      0   0   1 0.6369134    
## 26   0.14286      0   0   1  0.0952381      0   0   1 0.7521784    
## 27   0.00000      0   0   0  0.1095238      0   0   1 8.745e-07 ***
## 28   0.28571      0   0   1  0.2476190      0   0   1 0.8448380    
## 29   0.28571      0   0   1  0.1142857      0   0   1 0.3906779    
## 30   0.00000      0   0   0  0.0142857      0   0   1 0.0832610 .  
## 31   0.00000      0   0   0  0.0476190      0   0   1 0.0014252 ** 
## 32   0.14286      0   0   1  0.0666667      0   0   1 0.6149535    
## 33   0.00000      0   0   0  0.0333333      0   0   1 0.0078449 ** 
## 34   0.00000      0   0   0  0.0047619      0   0   1 0.3184669    
## 35   0.00000      0   0   0  0.0000000      0   0   0 1.0000000    
## 36   0.00000      0   0   0  0.0000000      0   0   0 1.0000000    
## 37   0.00000      0   0   0  0.0000000      0   0   0 1.0000000    
## 38   0.00000      0   0   0  0.0047619      0   0   1 0.3184669    
## 39   0.00000      0   0   0  0.1047619      0   0   1 1.559e-06 ***
## 40   0.14286      0   0   1  0.0285714      0   0   1 0.4552264    
## 41   0.00000      0   0   0  0.0047619      0   0   1 0.3184669    
## 42   0.00000      0   0   0  0.0000000      0   0   0 1.0000000    
## 43   0.00000      0   0   0  0.0000000      0   0   0 1.0000000    
## 44   0.28571      0   0   1  0.0809524      0   0   1 0.3108619    
## 45   0.28571      0   0   1  0.1571429      0   0   1 0.5146501    
## 46   0.00000      0   0   0  0.1238095      0   0   1 1.524e-07 ***
## 47   0.14286      0   0   1  0.1047619      0   0   1 0.8004085    
## 48   0.00000      0   0   0  0.0428571      0   0   1 0.0025105 ** 
## 49   0.14286      0   0   1  0.0714286      0   0   1 0.6369134    
## 50   0.14286      0   0   1  0.1714286      0   0   1 0.8501175    
## 51   0.14286      0   0   1  0.1904762      0   0   1 0.7536980    
## 52   0.28571      0   0   2  0.1047619      0   0   2 0.5508579    
## 53   0.00000      0   0   0  0.0809524      0   0   1 2.724e-05 ***
## 54   0.00000      0   0   0  0.0095238      0   0   1 0.1577941    
## 55   0.00000      0   0   0  0.0857143      0   0   1 1.542e-05 ***
## 56   0.14286      0   0   1  0.0571429      0   0   1 0.5722807    
## 57   0.14286      0   0   1  0.0047619      0   0   1 0.3712110    
## 58   0.00000      0   0   0  0.0000000      0   0   0 1.0000000

## 
## Note: Small p-values indicate a significant difference 
##       between simulations and observed network(s).

## Geodesic distances

##     obs: mean median  min   max  sim: mean median  min   max   Pr(>z)    
## 1      537.43    548  284   728 5.1929e+02    522  216   808 0.812902    
## 2     5018.86   4948 1252  7936 4.1557e+03   4342  608  8254 0.462454    
## 3     6750.86   5620 1542 11706 5.8890e+03   5059  980 12830 0.640328    
## 4     4814.00   5046 1550  7518 3.9995e+03   3968 1032  7100 0.423672    
## 5      471.71    536  116   812 8.7863e+02    974   68  1972 0.011941 *  
## 6       82.00     14    0   278 2.5990e+02    255    0   836 0.005932 ** 
## 7       10.00      0    0    38 8.0190e+01     30    0   472 6.11e-08 ***
## 8        0.00      0    0     0 1.9952e+01      2    0   272 2.33e-10 ***
## 9        0.00      0    0     0 5.2286e+00      0    0   152 2.81e-05 ***
## 10       0.00      0    0     0 1.3619e+00      0    0    88 0.009096 ** 
## 11       0.00      0    0     0 2.9524e-01      0    0    32 0.074865 .  
## 12       0.00      0    0     0 6.6667e-02      0    0     8 0.126931    
## 13       0.00      0    0     0 9.5238e-03      0    0     2 0.318467    
## 14       0.00      0    0     0 0.0000e+00      0    0     0 1.000000    
## Inf  14535.14  15184 4164 26670 1.6411e+04  17087 3828 28052 0.631090

## 
## Note: Small p-values indicate a significant difference 
##       between simulations and observed network(s).
```

Plot GOF statistics for `m0`

```r
plot(gof0)
```

![](data/amj_run_TERGM_tutorial_1_files/figure-html/gof_0_plot-1.png)<!-- -->

Now compare the GOF for `m1`

```r
gof1 <- gof(fit1, statistics=c(dsp, esp, deg, geodesic), nsim=30)
```

```
## 
## Starting GOF assessment on a single computing core....
## 
## Initial dimensions of the network and covariates:

##              t=2 t=3 t=4 t=5 t=6 t=7 t=8
## nets (row)   180 180 180 180 180 180 180
## nets (col)   180 180 180 180 180 180 180
## memory (row) 180 180 180 180 180 180 180
## memory (col) 180 180 180 180 180 180 180

## 
## All networks are conformable.
## 
## Dimensions of the network and covariates after adjustment:

##              t=2 t=3 t=4 t=5 t=6 t=7 t=8
## nets (row)   180 180 180 180 180 180 180
## nets (col)   180 180 180 180 180 180 180
## memory (row) 180 180 180 180 180 180 180
## memory (col) 180 180 180 180 180 180 180

## 
## No 'target' network(s) provided. Using networks on the left-hand side of the model formula as observed networks.
## Simulating 30 networks from the following formula:
##  nets[[1]] ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + nodematch("ipo_status", diff = F) + nodematch("state_code", diff = F) + nodecov("age") + absdiff("age") + nodecov("cent_deg") + edgecov(memory[[1]]) + nodecov("genidx_multilevel") + nodecov("cent_pow_n0_4") + absdiff("cent_pow_n0_4") + cycle(3) + cycle(4)
## Simulating 30 networks from the following formula:
##  nets[[2]] ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + nodematch("ipo_status", diff = F) + nodematch("state_code", diff = F) + nodecov("age") + absdiff("age") + nodecov("cent_deg") + edgecov(memory[[2]]) + nodecov("genidx_multilevel") + nodecov("cent_pow_n0_4") + absdiff("cent_pow_n0_4") + cycle(3) + cycle(4)
## Simulating 30 networks from the following formula:
##  nets[[3]] ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + nodematch("ipo_status", diff = F) + nodematch("state_code", diff = F) + nodecov("age") + absdiff("age") + nodecov("cent_deg") + edgecov(memory[[3]]) + nodecov("genidx_multilevel") + nodecov("cent_pow_n0_4") + absdiff("cent_pow_n0_4") + cycle(3) + cycle(4)
## Simulating 30 networks from the following formula:
##  nets[[4]] ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + nodematch("ipo_status", diff = F) + nodematch("state_code", diff = F) + nodecov("age") + absdiff("age") + nodecov("cent_deg") + edgecov(memory[[4]]) + nodecov("genidx_multilevel") + nodecov("cent_pow_n0_4") + absdiff("cent_pow_n0_4") + cycle(3) + cycle(4)
## Simulating 30 networks from the following formula:
##  nets[[5]] ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + nodematch("ipo_status", diff = F) + nodematch("state_code", diff = F) + nodecov("age") + absdiff("age") + nodecov("cent_deg") + edgecov(memory[[5]]) + nodecov("genidx_multilevel") + nodecov("cent_pow_n0_4") + absdiff("cent_pow_n0_4") + cycle(3) + cycle(4)
## Simulating 30 networks from the following formula:
##  nets[[6]] ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + nodematch("ipo_status", diff = F) + nodematch("state_code", diff = F) + nodecov("age") + absdiff("age") + nodecov("cent_deg") + edgecov(memory[[6]]) + nodecov("genidx_multilevel") + nodecov("cent_pow_n0_4") + absdiff("cent_pow_n0_4") + cycle(3) + cycle(4)
## Simulating 30 networks from the following formula:
##  nets[[7]] ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + nodematch("ipo_status", diff = F) + nodematch("state_code", diff = F) + nodecov("age") + absdiff("age") + nodecov("cent_deg") + edgecov(memory[[7]]) + nodecov("genidx_multilevel") + nodecov("cent_pow_n0_4") + absdiff("cent_pow_n0_4") + cycle(3) + cycle(4)
## 7 networks from which simulations are drawn were provided.
## Processing statistic: Dyad-wise shared partners
## Processing statistic: Edge-wise shared partners
## Processing statistic: Degree
## Processing statistic: Geodesic distances
```

```r
print(gof1)
```

```
## Dyad-wise shared partners

##     obs: mean median   min   max  sim: mean median   min   max    Pr(>z)
## 0  2.6815e+04  26890 23734 30776 2.6475e+04  25227 21022 31432 0.7812126
## 1  4.1563e+03   4294   972  6534 4.3773e+03   5625   600  8458 0.8167212
## 2  9.1171e+02    698   366  1456 9.8446e+02    991   130  2014 0.7169791
## 3  1.9686e+02    198    72   298 2.2584e+02    230    32   492 0.4910735
## 4  8.2857e+01     90    26   130 8.6505e+01     84     4   200 0.8287811
## 5  3.1143e+01     28     4    58 3.4676e+01     42     0    82 0.6901652
## 6  9.1429e+00      8     0    16 1.3457e+01     12     0    38 0.1716471
## 7  6.2857e+00      6     0    14 8.0476e+00      8     0    28 0.4362546
## 8  2.5714e+00      2     0     6 3.0857e+00      2     0    16 0.6471140
## 9  5.7143e-01      0     0     4 2.5429e+00      2     0    12 0.0125747
## 10 3.1429e+00      2     0     8 2.1810e+00      2     0    10 0.4900297
## 11 8.5714e-01      0     0     2 1.2381e+00      0     0    10 0.4002964
## 12 2.8571e-01      0     0     2 6.2857e-01      0     0     8 0.2891486
## 13 0.0000e+00      0     0     0 4.1905e-01      0     0     4 4.167e-10
## 14 8.5714e-01      0     0     2 3.5238e-01      0     0     6 0.2613795
## 15 0.0000e+00      0     0     0 3.1429e-01      0     0     4 2.010e-08
## 16 0.0000e+00      0     0     0 2.8571e-01      0     0     4 4.364e-08
## 17 0.0000e+00      0     0     0 2.6667e-01      0     0     4 1.409e-07
## 18 2.8571e-01      0     0     2 2.6667e-01      0     0     2 0.9495807
## 19 2.8571e-01      0     0     2 2.9524e-01      0     0     2 0.9747995
## 20 0.0000e+00      0     0     0 3.2381e-01      0     0     2 1.291e-09
## 21 2.8571e-01      0     0     2 2.7619e-01      0     0     2 0.9747830
## 22 5.7143e-01      0     0     2 1.6190e-01      0     0     2 0.3108619
## 23 2.8571e-01      0     0     2 1.9048e-01      0     0     4 0.7523925
## 24 2.8571e-01      0     0     2 2.5714e-01      0     0     4 0.9246306
## 25 8.5714e-01      0     0     2 4.0000e-01      0     0     4 0.3039967
## 26 0.0000e+00      0     0     0 3.6190e-01      0     0     4 1.364e-08
## 27 0.0000e+00      0     0     0 3.3333e-01      0     0     4 3.519e-08
## 28 0.0000e+00      0     0     0 1.0476e-01      0     0     2 0.0008101
## 29 0.0000e+00      0     0     0 4.7619e-02      0     0     2 0.0249929
## 30 0.0000e+00      0     0     0 3.8095e-02      0     0     2 0.0452374
## 31 0.0000e+00      0     0     0 9.5238e-03      0     0     2 0.3184669
## 32 0.0000e+00      0     0     0 9.5238e-03      0     0     2 0.3184669
## 33 0.0000e+00      0     0     0 0.0000e+00      0     0     0 1.0000000
##       
## 0     
## 1     
## 2     
## 3     
## 4     
## 5     
## 6     
## 7     
## 8     
## 9  *  
## 10    
## 11    
## 12    
## 13 ***
## 14    
## 15 ***
## 16 ***
## 17 ***
## 18    
## 19    
## 20 ***
## 21    
## 22    
## 23    
## 24    
## 25    
## 26 ***
## 27 ***
## 28 ***
## 29 *  
## 30 *  
## 31    
## 32    
## 33

## 
## Note: Small p-values indicate a significant difference 
##       between simulations and observed network(s).

## Edge-wise shared partners

##    obs: mean median min max  sim: mean median min max    Pr(>z)    
## 0  150.85714    166  92 194 1.3926e+02    138  70 194 0.4542258    
## 1  167.14286    156 102 234 1.6317e+02    169  60 284 0.8604042    
## 2   84.57143     78  42 126 9.1562e+01     97  14 176 0.6043855    
## 3   62.85714     74  30  92 5.8486e+01     50  18 120 0.6784580    
## 4   32.85714     36  14  54 3.5238e+01     35   0  78 0.6899464    
## 5   18.57143     18   2  34 1.9067e+01     20   0  48 0.9283142    
## 6    7.14286      8   0  12 9.8095e+00     10   0  28 0.2493395    
## 7    4.57143      4   0  10 6.5714e+00      6   0  24 0.2557805    
## 8    1.71429      2   0   4 2.2571e+00      2   0  10 0.4681195    
## 9    0.57143      0   0   4 2.2000e+00      2   0  10 0.0288333 *  
## 10   3.14286      2   0   8 2.0667e+00      1   0  10 0.4420056    
## 11   0.57143      0   0   2 1.2381e+00      0   0  10 0.1300864    
## 12   0.00000      0   0   0 6.2857e-01      0   0   8 6.361e-11 ***
## 13   0.00000      0   0   0 4.1905e-01      0   0   4 4.167e-10 ***
## 14   0.00000      0   0   0 3.4286e-01      0   0   6 2.860e-07 ***
## 15   0.00000      0   0   0 2.6667e-01      0   0   4 3.660e-07 ***
## 16   0.00000      0   0   0 1.8095e-01      0   0   4 2.445e-05 ***
## 17   0.00000      0   0   0 1.3333e-01      0   0   2 0.0001490 ***
## 18   0.28571      0   0   2 1.1429e-01      0   0   2 0.5722807    
## 19   0.28571      0   0   2 1.4286e-01      0   0   2 0.6369134    
## 20   0.00000      0   0   0 2.3810e-01      0   0   2 2.735e-07 ***
## 21   0.28571      0   0   2 2.1905e-01      0   0   2 0.8248759    
## 22   0.57143      0   0   2 1.4286e-01      0   0   2 0.2906634    
## 23   0.28571      0   0   2 1.5238e-01      0   0   4 0.6595427    
## 24   0.28571      0   0   2 2.3810e-01      0   0   4 0.8746005    
## 25   0.85714      0   0   2 3.9048e-01      0   0   4 0.2948466    
## 26   0.00000      0   0   0 3.6190e-01      0   0   4 1.364e-08 ***
## 27   0.00000      0   0   0 3.3333e-01      0   0   4 3.519e-08 ***
## 28   0.00000      0   0   0 1.0476e-01      0   0   2 0.0008101 ***
## 29   0.00000      0   0   0 4.7619e-02      0   0   2 0.0249929 *  
## 30   0.00000      0   0   0 3.8095e-02      0   0   2 0.0452374 *  
## 31   0.00000      0   0   0 9.5238e-03      0   0   2 0.3184669    
## 32   0.00000      0   0   0 9.5238e-03      0   0   2 0.3184669    
## 33   0.00000      0   0   0 0.0000e+00      0   0   0 1.0000000

## 
## Note: Small p-values indicate a significant difference 
##       between simulations and observed network(s).

## Degree

##    obs: mean median min max  sim: mean median min max    Pr(>z)    
## 0   50.71429     46  12 105 55.2714286   48.0   9 118 0.7702925    
## 1   45.00000     45  22  61 38.9476190   38.0  22  58 0.3550056    
## 2   21.28571     24  14  26 24.7238095   26.0   9  41 0.1469390    
## 3   18.28571     16  11  25 16.6809524   17.0   2  34 0.5047365    
## 4    9.42857      9   7  12 10.0619048   10.5   1  21 0.4557195    
## 5   11.42857     12   6  17  9.5666667    9.0   1  21 0.2719355    
## 6    8.57143     10   4  12  7.3523810    7.0   1  16 0.3205023    
## 7    2.00000      2   1   3  4.3761905    4.0   0  10 6.348e-05 ***
## 8    3.71429      4   3   5  2.7523810    3.0   0   8 0.0140353 *  
## 9    1.42857      2   0   2  1.8952381    2.0   0   4 0.1741906    
## 10   1.00000      1   0   2  1.3380952    1.0   0   5 0.3231522    
## 11   1.42857      1   1   2  0.9714286    1.0   0   5 0.0665010 .  
## 12   0.42857      0   0   1  0.7238095    1.0   0   3 0.2010782    
## 13   0.42857      0   0   1  0.5285714    0.0   0   4 0.6467909    
## 14   0.28571      0   0   1  0.2809524    0.0   0   2 0.9805549    
## 15   0.28571      0   0   1  0.3476190    0.0   0   2 0.7529167    
## 16   0.00000      0   0   0  0.3714286    0.0   0   2 < 2.2e-16 ***
## 17   0.00000      0   0   0  0.1809524    0.0   0   1 1.108e-10 ***
## 18   0.28571      0   0   1  0.1428571    0.0   0   2 0.4710908    
## 19   0.28571      0   0   1  0.1809524    0.0   0   2 0.5937300    
## 20   0.14286      0   0   1  0.1190476    0.0   0   2 0.8744872    
## 21   0.14286      0   0   1  0.2000000    0.0   0   2 0.7076634    
## 22   0.57143      1   0   1  0.1476190    0.0   0   2 0.0811424 .  
## 23   0.00000      0   0   0  0.1857143    0.0   0   2 5.591e-10 ***
## 24   0.00000      0   0   0  0.1285714    0.0   0   1 8.473e-08 ***
## 25   0.14286      0   0   1  0.1047619    0.0   0   1 0.8004085    
## 26   0.14286      0   0   1  0.0619048    0.0   0   1 0.5934037    
## 27   0.00000      0   0   0  0.1047619    0.0   0   1 1.559e-06 ***
## 28   0.28571      0   0   1  0.1714286    0.0   0   2 0.5613693    
## 29   0.28571      0   0   1  0.1285714    0.0   0   1 0.4293359    
## 30   0.00000      0   0   0  0.1619048    0.0   0   1 1.291e-09 ***
## 31   0.00000      0   0   0  0.0952381    0.0   0   1 4.921e-06 ***
## 32   0.14286      0   0   1  0.0666667    0.0   0   1 0.6149535    
## 33   0.00000      0   0   0  0.0285714    0.0   0   1 0.0139540 *  
## 34   0.00000      0   0   0  0.0095238    0.0   0   1 0.1577941    
## 35   0.00000      0   0   0  0.0142857    0.0   0   1 0.0832610 .  
## 36   0.00000      0   0   0  0.0333333    0.0   0   1 0.0078449 ** 
## 37   0.00000      0   0   0  0.0428571    0.0   0   1 0.0025105 ** 
## 38   0.00000      0   0   0  0.0190476    0.0   0   1 0.0452374 *  
## 39   0.00000      0   0   0  0.0380952    0.0   0   1 0.0044309 ** 
## 40   0.14286      0   0   1  0.0000000    0.0   0   0 0.3559177    
## 41   0.00000      0   0   0  0.0095238    0.0   0   1 0.1577941    
## 42   0.00000      0   0   0  0.0047619    0.0   0   1 0.3184669    
## 43   0.00000      0   0   0  0.0000000    0.0   0   0 1.0000000    
## 44   0.28571      0   0   1  0.0047619    0.0   0   1 0.1785568    
## 45   0.28571      0   0   1  0.0238095    0.0   0   1 0.2057271    
## 46   0.00000      0   0   0  0.0428571    0.0   0   1 0.0025105 ** 
## 47   0.14286      0   0   1  0.0333333    0.0   0   1 0.4735323    
## 48   0.00000      0   0   0  0.0380952    0.0   0   1 0.0044309 ** 
## 49   0.14286      0   0   1  0.0333333    0.0   0   1 0.4735323    
## 50   0.14286      0   0   1  0.0428571    0.0   0   1 0.5116132    
## 51   0.14286      0   0   1  0.0809524    0.0   0   1 0.6819897    
## 52   0.28571      0   0   2  0.0714286    0.0   0   1 0.4822211    
## 53   0.00000      0   0   0  0.1190476    0.0   0   2 2.018e-06 ***
## 54   0.00000      0   0   0  0.1047619    0.0   0   2 4.473e-06 ***
## 55   0.00000      0   0   0  0.1142857    0.0   0   2 7.778e-06 ***
## 56   0.14286      0   0   1  0.0904762    0.0   0   1 0.7284695    
## 57   0.14286      0   0   1  0.0619048    0.0   0   2 0.5937159    
## 58   0.00000      0   0   0  0.0523810    0.0   0   1 0.0008101 ***
## 59   0.00000      0   0   0  0.0380952    0.0   0   1 0.0044309 ** 
## 60   0.00000      0   0   0  0.0238095    0.0   0   1 0.0249929 *  
## 61   0.00000      0   0   0  0.0761905    0.0   0   1 4.804e-05 ***
## 62   0.00000      0   0   0  0.0714286    0.0   0   2 0.0002294 ***
## 63   0.00000      0   0   0  0.0666667    0.0   0   1 0.0001490 ***
## 64   0.00000      0   0   0  0.0761905    0.0   0   2 0.0001314 ***
## 65   0.00000      0   0   0  0.0380952    0.0   0   1 0.0044309 ** 
## 66   0.00000      0   0   0  0.0333333    0.0   0   1 0.0078449 ** 
## 67   0.00000      0   0   0  0.0095238    0.0   0   1 0.1577941    
## 68   0.00000      0   0   0  0.0047619    0.0   0   1 0.3184669    
## 69   0.00000      0   0   0  0.0190476    0.0   0   1 0.0452374 *  
## 70   0.00000      0   0   0  0.0000000    0.0   0   0 1.0000000    
## 71   0.00000      0   0   0  0.0047619    0.0   0   1 0.3184669    
## 72   0.00000      0   0   0  0.0000000    0.0   0   0 1.0000000    
## 73   0.00000      0   0   0  0.0095238    0.0   0   1 0.1577941    
## 74   0.00000      0   0   0  0.0095238    0.0   0   1 0.1577941    
## 75   0.00000      0   0   0  0.0047619    0.0   0   1 0.3184669    
## 76   0.00000      0   0   0  0.0095238    0.0   0   1 0.1577941    
## 77   0.00000      0   0   0  0.0095238    0.0   0   1 0.1577941    
## 78   0.00000      0   0   0  0.0095238    0.0   0   1 0.1577941    
## 79   0.00000      0   0   0  0.0000000    0.0   0   0 1.0000000

## 
## Note: Small p-values indicate a significant difference 
##       between simulations and observed network(s).

## Geodesic distances

##     obs: mean median  min   max  sim: mean median  min   max   Pr(>z)   
## 1      537.43    548  284   728 5.3544e+02    573  204   808 0.979255   
## 2     5018.86   4948 1252  7936 5.3485e+03   6567  658 10568 0.776198   
## 3     6750.86   5620 1542 11706 6.4890e+03   6130  952 13492 0.886231   
## 4     4814.00   5046 1550  7518 3.2026e+03   3027  944  6400 0.138854   
## 5      471.71    536  116   812 5.4467e+02    558   20  1294 0.559407   
## 6       82.00     14    0   278 1.1897e+02     77    0   450 0.442898   
## 7       10.00      0    0    38 1.9933e+01      2    0   256 0.191117   
## 8        0.00      0    0     0 2.4381e+00      0    0   120 0.001334 **
## 9        0.00      0    0     0 2.3810e-01      0    0     8 0.005793 **
## 10       0.00      0    0     0 1.9048e-02      0    0     2 0.157794   
## 11       0.00      0    0     0 0.0000e+00      0    0     0 1.000000   
## Inf  14535.14  15184 4164 26670 1.5958e+04  15448 3150 28678 0.714688

## 
## Note: Small p-values indicate a significant difference 
##       between simulations and observed network(s).
```

Plot GOF statistics for `m1`

```r
plot(gof1)
```

![](data/amj_run_TERGM_tutorial_1_files/figure-html/gof_1_plot-1.png)<!-- -->



# Part 6: Estimation Algorithm

> [Back to Contents](#contents  "Back")

Compare the estamation algorithm for the model `m4`.  You may repeat the same steps as needed to compare other models.

Compute the first model `m4` again using MCMCMLE (instead of bootstrapped PMLE) and save to disk as an RDS (serialized) file:

```r
## set pseudorandom number generator seed for reproducibility
set.seed(1111)
## estimate the TERGM with bootstrapped PMLE
fit0m <- mtergm(m0, ctrl=control.ergm(seed = 1111))
```

```
## 
## Initial dimensions of the network and covariates:
##              t=2 t=3 t=4 t=5 t=6 t=7 t=8
## nets (row)   180 180 180 180 180 180 180
## nets (col)   180 180 180 180 180 180 180
## memory (row) 180 180 180 180 180 180 180
## memory (col) 180 180 180 180 180 180 180
## 
## All networks are conformable.
## 
## Dimensions of the network and covariates after adjustment:
##              t=2 t=3 t=4 t=5 t=6 t=7 t=8
## nets (row)   180 180 180 180 180 180 180
## nets (col)   180 180 180 180 180 180 180
## memory (row) 180 180 180 180 180 180 180
## memory (col) 180 180 180 180 180 180 180

## Estimating...
## Starting maximum likelihood estimation via MCMLE:
## Iteration 1 of at most 20:
## Optimizing with step length 0.380849430581003.
## The log-likelihood improved by 2.465.
## Iteration 2 of at most 20:
## Optimizing with step length 0.659003375233084.
## The log-likelihood improved by 2.668.
## Iteration 3 of at most 20:
## Optimizing with step length 1.
## The log-likelihood improved by 1.182.
## Step length converged once. Increasing MCMC sample size.
## Iteration 4 of at most 20:
## Optimizing with step length 1.
## The log-likelihood improved by 1.42.
## Step length converged twice. Stopping.
## Evaluating log-likelihood at the estimate. Using 20 bridges: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 .
## This model was fit using MCMC.  To examine model diagnostics and check for degeneracy, use the mcmc.diagnostics() function.
## Done.
```

```r
## SAVE SERIALIZED DATA
fit0m_file <- file.path(data_dir,sprintf('fit_%s_pd%s_%s.rds', firm_i, nPeriods, 'm0m'))
saveRDS(fit0m, file=fit0m_file)
```

Compare PMLE and MCMCMLE results

```r
## Cache model fits list
fits <- list(PMLE=fit0, MCMCMLE=fit0m)

## Echo model comparison table to screen
screenreg(fits, digits = 3)
```

```
## 
## ======================================================
##                       PMLE              MCMCMLE       
## ------------------------------------------------------
## edges                    -0.180             -0.477    
##                       [-1.749;  1.067]      (0.382)   
## gwesp.fixed.0             0.813 *            0.964 ***
##                       [ 0.658;  0.908]      (0.076)   
## gwdegree                 -1.145 *           -0.741 ***
##                       [-1.910; -0.527]      (0.144)   
## nodematch.ipo_status      0.132              0.175    
##                       [-0.259;  1.070]      (0.358)   
## nodematch.state_code     -0.593 *           -0.624 ***
##                       [-0.770; -0.387]      (0.169)   
## nodecov.age              -0.139 *           -0.128 ***
##                       [-0.199; -0.065]      (0.011)   
## absdiff.age               0.157 *            0.145 ***
##                       [ 0.085;  0.216]      (0.012)   
## edgecov.memory[[i]]       5.229 *                     
##                       [ 4.908;  5.769]                
## edgecov.memory                               5.047 ***
##                                             (0.144)   
## ------------------------------------------------------
## Num. obs.             13787             225540        
## ======================================================
## *** p < 0.001, ** p < 0.01, * p < 0.05 (or 0 outside the confidence interval).
```

```r
## SAVE FORMATTED REGRESSION TABLE
compare_file <- file.path(data_dir,sprintf('%s_tergm_results_pd%s_R%s_%s.html', firm_i, nPeriods, R, 'm0PMLE-m0MCMCMLE'))
htmlreg(fits, digits = 3, file=compare_file)
```

And compare the PMLE and MCMCMLE confidence intervals directly

```r
## Cache model fits list
fits <- list(PMLE=fit0, MCMCMLE=fit0m)

## Echo model comparison table to screen
screenreg(fits, digits = 3, ci.force = T, ci.force.level = .95)
```

```
## 
## ========================================================
##                       PMLE              MCMCMLE         
## --------------------------------------------------------
## edges                    -0.180             -0.477      
##                       [-1.749;  1.067]  [-1.225;  0.271]
## gwesp.fixed.0             0.813 *            0.964 *    
##                       [ 0.658;  0.908]  [ 0.814;  1.113]
## gwdegree                 -1.145 *           -0.741 *    
##                       [-1.910; -0.527]  [-1.023; -0.459]
## nodematch.ipo_status      0.132              0.175      
##                       [-0.259;  1.070]  [-0.527;  0.878]
## nodematch.state_code     -0.593 *           -0.624 *    
##                       [-0.770; -0.387]  [-0.956; -0.292]
## nodecov.age              -0.139 *           -0.128 *    
##                       [-0.199; -0.065]  [-0.150; -0.107]
## absdiff.age               0.157 *            0.145 *    
##                       [ 0.085;  0.216]  [ 0.122;  0.168]
## edgecov.memory[[i]]       5.229 *                       
##                       [ 4.908;  5.769]                  
## edgecov.memory                               5.047 *    
##                                         [ 4.765;  5.330]
## --------------------------------------------------------
## Num. obs.             13787             225540          
## ========================================================
## * 0 outside the confidence interval
```

```r
## SAVE FORMATTED REGRESSION TABLE
compare_file <- file.path(data_dir,sprintf('%s_tergm_results_pd%s_R%s_%s.html', firm_i, nPeriods, R, 'm0PMLE-m0MCMCMLE_ci'))
htmlreg(fits, digits = 3, file=compare_file, ci.force = T, ci.force.level = .05)
```

Finally, check the MCMCMLE diagnostics

```r
mcmc.diagnostics(fit0m@ergm)
```

```
## Sample statistics summary:
## 
## Iterations = 16384:4209664
## Thinning interval = 1024 
## Number of chains = 1 
## Sample size per chain = 4096 
## 
## 1. Empirical mean and standard deviation for each variable,
##    plus standard error of the mean:
## 
##                           Mean      SD Naive SE Time-series SE
## edges                   3.0723  19.700   0.3078         1.2698
## gwesp.fixed.0           3.5552  19.959   0.3119         1.8134
## gwdegree               -0.3508  11.232   0.1755         0.5734
## nodematch.ipo_status    1.2542  19.092   0.2983         1.3942
## nodematch.state_code   -3.6062   6.416   0.1002         0.3763
## nodecov.age           -35.6899 384.359   6.0056        27.4444
## absdiff.age          -134.2041 300.434   4.6943        19.9207
## edgecov.memory         10.8696  19.003   0.2969         1.1498
## 
## 2. Quantiles for each variable:
## 
##                         2.5%  25%  50% 75% 97.5%
## edges                 -36.00  -10    3  17  42.0
## gwesp.fixed.0         -39.00   -9    4  17  41.0
## gwdegree              -21.62   -8   -1   7  22.0
## nodematch.ipo_status  -36.62  -11    1  14  38.0
## nodematch.state_code  -16.00   -8   -4   1  10.0
## nodecov.age          -784.00 -307  -28 231 685.6
## absdiff.age          -702.00 -343 -132  73 453.6
## edgecov.memory        -27.00   -2   12  24  47.0
## 
## 
## Sample statistics cross-correlations:
##                           edges gwesp.fixed.0   gwdegree
## edges                 1.0000000     0.6314057  0.6581411
## gwesp.fixed.0         0.6314057     1.0000000  0.1204709
## gwdegree              0.6581411     0.1204709  1.0000000
## nodematch.ipo_status  0.9886884     0.6256350  0.6462060
## nodematch.state_code  0.3580458     0.2441033  0.2293714
## nodecov.age           0.7275301     0.4828927  0.3857011
## absdiff.age           0.5918817     0.3200693  0.3720640
## edgecov.memory       -0.8724774    -0.5287060 -0.5860218
##                      nodematch.ipo_status nodematch.state_code nodecov.age
## edges                           0.9886884            0.3580458   0.7275301
## gwesp.fixed.0                   0.6256350            0.2441033   0.4828927
## gwdegree                        0.6462060            0.2293714   0.3857011
## nodematch.ipo_status            1.0000000            0.3644899   0.7101860
## nodematch.state_code            0.3644899            1.0000000   0.2704267
## nodecov.age                     0.7101860            0.2704267   1.0000000
## absdiff.age                     0.5755261            0.2028273   0.9265124
## edgecov.memory                 -0.8753191           -0.3510019  -0.5997116
##                      absdiff.age edgecov.memory
## edges                  0.5918817     -0.8724774
## gwesp.fixed.0          0.3200693     -0.5287060
## gwdegree               0.3720640     -0.5860218
## nodematch.ipo_status   0.5755261     -0.8753191
## nodematch.state_code   0.2028273     -0.3510019
## nodecov.age            0.9265124     -0.5997116
## absdiff.age            1.0000000     -0.5667370
## edgecov.memory        -0.5667370      1.0000000
## 
## Sample statistics auto-correlation:
## Chain 1 
##              edges gwesp.fixed.0  gwdegree nodematch.ipo_status
## Lag 0    1.0000000     1.0000000 1.0000000            1.0000000
## Lag 1024 0.8838772     0.9425326 0.8096165            0.8798159
## Lag 2048 0.7864046     0.8876284 0.6656872            0.7803422
## Lag 3072 0.7030145     0.8392670 0.5566241            0.6960648
## Lag 4096 0.6326490     0.7922254 0.4692374            0.6253648
## Lag 5120 0.5733575     0.7489398 0.3961800            0.5646937
##          nodematch.state_code nodecov.age absdiff.age edgecov.memory
## Lag 0               1.0000000   1.0000000   1.0000000      1.0000000
## Lag 1024            0.8347718   0.8881187   0.8874959      0.8749346
## Lag 2048            0.7022984   0.7968271   0.7951130      0.7698871
## Lag 3072            0.5937905   0.7178424   0.7125952      0.6797544
## Lag 4096            0.4979769   0.6498921   0.6418456      0.6035623
## Lag 5120            0.4202840   0.5955045   0.5839419      0.5391573
## 
## Sample statistics burn-in diagnostic (Geweke):
## Chain 1 
## 
## Fraction in 1st window = 0.1
## Fraction in 2nd window = 0.5 
## 
##                edges        gwesp.fixed.0             gwdegree 
##              -2.4150              -1.6614              -1.7590 
## nodematch.ipo_status nodematch.state_code          nodecov.age 
##              -1.8425               1.1636              -1.2394 
##          absdiff.age       edgecov.memory 
##              -1.3510              -0.3909 
## 
## Individual P-values (lower = worse):
##                edges        gwesp.fixed.0             gwdegree 
##           0.01573314           0.09663062           0.07857317 
## nodematch.ipo_status nodematch.state_code          nodecov.age 
##           0.06539507           0.24456780           0.21520396 
##          absdiff.age       edgecov.memory 
##           0.17668597           0.69585895 
## Joint P-value (lower = worse):  5.638501e-08 .

## Warning in formals(fun): argument is not a function
```

![](data/amj_run_TERGM_tutorial_1_files/figure-html/mcmc_diag-1.png)<!-- -->![](data/amj_run_TERGM_tutorial_1_files/figure-html/mcmc_diag-2.png)<!-- -->![](data/amj_run_TERGM_tutorial_1_files/figure-html/mcmc_diag-3.png)<!-- -->

```
## 
## MCMC diagnostics shown here are from the last round of simulation, prior to computation of final parameter estimates. Because the final estimates are refinements of those used for this simulation run, these diagnostics may understate model performance. To directly assess the performance of the final model on in-model statistics, please use the GOF command: gof(ergmFitObject, GOF=~model).
```

That concludes the main analyses.

