02a-flowEMMi for Z-project
================
Compiled at 2023-08-22 00:00:55 UTC

``` r
library("conflicted")
library(purrr)
library(dplyr)
library(mvtnorm)
library(flowEMMi)
library(flowCore)
library(flowWorkspace)
library(ggcyto)
library(ggforce)
library(tidyverse)
library(RColorBrewer)
library(knitr)
```

## Z-project

### Import preprocessed data

``` r
DAPI <- readRDS("~/Desktop/MScThesis/workflow/data/new_DAPI.rds")

FDA_PI <- readRDS("~/Desktop/MScThesis/workflow/data/new_FDA_PI.rds")

gating_DAPI <- readRDS("~/Desktop/z_gating_list.rds")
```

### flowEMMI gating on DAPI

Based on the scatter plot in 01-data, we could tell that except for the
sample taken in the surrounding region, the cells of other region
samples are mostly concentrated in the range of 20000 to 60000.

Therefore, we will run the automated gating on surrounding region and
other regions respectively.

-   Setting range for other regions: 20000-50000

-   Setting range for Surrounding region: 1000-50000

``` r
set.seed(1)
location <- c("Inner_zone","Middle_zone","Outer_zone","Whole_colony","Surrounding")
gating_DAPI <- list()

for (i in 1:4){
  data_name <- paste0(location[i],"_DAPI.fcs")
  data <- DAPI[[i]]
  fdo <- mkFlowDataObject(data, xChannel="PMT.1", yChannel="PMT.9")
  gating <- flowEMMi( fdo=fdo, xMin=20000, xMax=50000, yMin=20000, yMax=50000
                      , initFraction=0.01
                      , finalFraction=1.0
                      , minClusters=5, maxClusters=15, clusterbracket=2
                      , numberOfInits=5
                      , verbose=TRUE
                      , parallel=FALSE
                      , convergenceEpsilon=0.01
                      , whenToRemoveOverlaps = 20
                      , mergeWhenCenter = FALSE
                      , mergeWhenTwoCenters = FALSE
                      , thresholdForDeletion = 0.2
                      , threshold = 0.9
                      , considerWeights=TRUE
                      , plot = FALSE
                      , alpha=0.9
                      , minMinor=500)
  gating_DAPI[[i]] <- gating$best
}

# DAPI surrounding
DAPI_sur <- DAPI[["Surrounding_DAPI.fcs"]] 
fdo_sur <- mkFlowDataObject(DAPI_sur, xChannel="PMT.1", yChannel="PMT.9")
gating_sur <- flowEMMi( fdo=fdo_sur, xMin=1000, xMax=50000, yMin=1000, yMax=50000
                      , initFraction=0.01
                      , finalFraction=1.0
                      , minClusters=5, maxClusters=15, clusterbracket=2
                      , numberOfInits=5
                      , verbose=TRUE
                      , parallel=FALSE
                      , convergenceEpsilon=0.01
                      , whenToRemoveOverlaps = 20
                      , mergeWhenCenter = FALSE
                      , mergeWhenTwoCenters = FALSE
                      , thresholdForDeletion = 0.2
                      , threshold = 0.9
                      , considerWeights=TRUE
                      , plot = FALSE
                      , alpha=0.9
                      , minMinor=500)
gating_DAPI[[5]] <- gating_sur$best
```

### Gating plots on DAPI

``` r
location <- c("Inner_zone","Middle_zone","Outer_zone","Whole_colony","Surrounding")
gating_DAPI_plot <- list()

for (i in 1:5){
  data_name <- paste0(location[i],"_DAPI.fcs")
  data <- DAPI[[i]]
  plots <- plotDensityAndEllipses(fcsData = data, ch1="PMT.1", ch2="PMT.9", alpha=0.9,
                            logScale = F, results = gating_DAPI[[i]],
                            title = data_name, plotRelevance = T,
                            ellipseDotSize = 0.5, axis_size=10, axisLabeling_size=10,
                            xlab = "Forward Scatter", ylab = "DAPI", font = "Arial")
  gating_DAPI_plot[[i]] <- plots$plot
}
```

    ## Warning in KernSmooth::bkde2D(x, bandwidth = bandwidth, gridsize = nbin, :
    ## Binning grid too coarse for current (small) bandwidth: consider increasing
    ## 'gridsize'

![](02a-flowEMMi-for-Z-project_files/figure-gfm/plot-DAPI-1.png)<!-- -->

    ## Warning in KernSmooth::bkde2D(x, bandwidth = bandwidth, gridsize = nbin, :
    ## Binning grid too coarse for current (small) bandwidth: consider increasing
    ## 'gridsize'

![](02a-flowEMMi-for-Z-project_files/figure-gfm/plot-DAPI-2.png)<!-- -->![](02a-flowEMMi-for-Z-project_files/figure-gfm/plot-DAPI-3.png)<!-- -->![](02a-flowEMMi-for-Z-project_files/figure-gfm/plot-DAPI-4.png)<!-- -->![](02a-flowEMMi-for-Z-project_files/figure-gfm/plot-DAPI-5.png)<!-- -->

In the above graphs, the ellipses are colored according to their cluster
probability, high probability in red, medium probability in white, low
probability in blue.

### Mahalanobis distance

``` r
maha_result <- list()

for (x in 1:5){
  data <- DAPI[[x]]@exprs[,c(11,27)]
  
  mu <- split(gating_DAPI[[x]]@mu,seq_len(ncol(gating_DAPI[[x]]@mu)))
  sigma <- gating_DAPI[[x]]@sigma
  
  n_cells <- nrow(data)
  n_cluster <- length(sigma)
  maha_data <- matrix(NA,nrow=n_cells,ncol=n_cluster)
  
  for (i in 1:n_cells){
    for (j in 1:n_cluster){
      maha_data[i,j] <- mahalanobis(data[i,],mu[[j]],sigma[[j]])
    }
  }
  maha_result[[x]] <- maha_data
}


# rename distance matrix
for (i in 1:5){
  maha1 <- maha_result[[i]]
  coln <- ncol(maha1)
  maha <- maha1[,2:coln] %>% as.data.frame()
  
  # set 50 as distance cutoff
  for (j in 1:nrow(maha)){
    rv <- maha[j,1:ncol(maha)-1]
    if(all(rv>50)) { maha$Min[j] <- NA}
    else {maha$Min[j] <- which.min(rv)}
  }

  maha_result[[i]] <- maha
  
  table <-table(maha[,ncol(maha)]) %>% 
    kable(caption = paste0(location[i],"_DAPI"),
          col.names = c("No. of Cluster","Number of Cells"))
  
  print(table)
}
```

    ## 
    ## 
    ## Table: Inner_zone_DAPI
    ## 
    ## |No. of Cluster | Number of Cells|
    ## |:--------------|---------------:|
    ## |1              |             166|
    ## |3              |             103|
    ## |4              |           42485|
    ## |5              |             848|
    ## 
    ## 
    ## Table: Middle_zone_DAPI
    ## 
    ## |No. of Cluster | Number of Cells|
    ## |:--------------|---------------:|
    ## |2              |           37248|
    ## |3              |             571|
    ## |4              |              13|
    ## |6              |             183|
    ## |7              |            2745|
    ## |8              |            9143|
    ## 
    ## 
    ## Table: Outer_zone_DAPI
    ## 
    ## |No. of Cluster | Number of Cells|
    ## |:--------------|---------------:|
    ## |1              |           35059|
    ## |2              |             239|
    ## |3              |            7729|
    ## 
    ## 
    ## Table: Whole_colony_DAPI
    ## 
    ## |No. of Cluster | Number of Cells|
    ## |:--------------|---------------:|
    ## |1              |            7467|
    ## |2              |           18573|
    ## |3              |               2|
    ## |4              |           16260|
    ## |5              |           19211|
    ## 
    ## 
    ## Table: Surrounding_DAPI
    ## 
    ## |No. of Cluster | Number of Cells|
    ## |:--------------|---------------:|
    ## |1              |             660|
    ## |2              |           35549|

``` r
# plot
DAPI_maha_plot <- list()

for (i in 1:5){
  data1 <- DAPI[[i]]@exprs[,c(11,27)] 
  data2 <- maha_result[[i]]
  data <- cbind(data1,data2$Min) %>% as.data.frame()
  colnames(data)[3] <- "Cluster"
  data$Cluster <- as.factor(data$Cluster)
  
  #ellipses generated by 
  plot1 <- ggplot(data,aes(x=PMT.1,y=PMT.9))+
    geom_point(aes(color=Cluster))+geom_mark_ellipse(aes(color=Cluster))+
    ggtitle(paste0(location[i],"_DAPI"),subtitle = "Auto-generated ellipses")
  
  print(plot1)
  
  #ellipses from flowEMMI result
  ellipses <- gating_DAPI[[i]]
  num_ellipse <- length(ellipses@sigma)
  
  p <- plot(NULL,type="n",xlim=c(0,70000),ylim=c(0,70000),xlab="PMT.1",ylab="PMT.9")
  
  for (j in 2:num_ellipse){
    mu <- ellipses@mu[,j]
    sigma <- ellipses@sigma[[j]]
    eli <- ellipse::ellipse(centre=mu,x=sigma,level=0.95,npoints=200) 
    p <- p+lines(eli,type="l",col=j)
  }
}
```

    ## Warning: Using the `size` aesthetic in this geom was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` in the `default_aes` field and elsewhere instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](02a-flowEMMi-for-Z-project_files/figure-gfm/maha-plot-1.png)<!-- -->![](02a-flowEMMi-for-Z-project_files/figure-gfm/maha-plot-2.png)<!-- -->![](02a-flowEMMi-for-Z-project_files/figure-gfm/maha-plot-3.png)<!-- -->![](02a-flowEMMi-for-Z-project_files/figure-gfm/maha-plot-4.png)<!-- -->![](02a-flowEMMi-for-Z-project_files/figure-gfm/maha-plot-5.png)<!-- -->![](02a-flowEMMi-for-Z-project_files/figure-gfm/maha-plot-6.png)<!-- -->![](02a-flowEMMi-for-Z-project_files/figure-gfm/maha-plot-7.png)<!-- -->![](02a-flowEMMi-for-Z-project_files/figure-gfm/maha-plot-8.png)<!-- -->![](02a-flowEMMi-for-Z-project_files/figure-gfm/maha-plot-9.png)<!-- -->![](02a-flowEMMi-for-Z-project_files/figure-gfm/maha-plot-10.png)<!-- -->

## Files written

These files have been written to the target directory,
`data/02a-flowEMMi for Z-project`:

``` r
projthis::proj_dir_info(path_target())
```

    ## # A tibble: 0 × 4
    ## # ℹ 4 variables: path <fs::path>, type <fct>, size <fs::bytes>,
    ## #   modification_time <dttm>