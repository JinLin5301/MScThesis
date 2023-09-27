02b-flowEMMI on FDA_PI from Z-project
================
Compiled at 2023-09-27 21:58:51 UTC

``` r
here::i_am(paste0(params$name, ".Rmd"), uuid = "ce529353-7643-4521-8dc9-35a25f812de1")
```

``` r
library("conflicted")
library(purrr)
library(dplyr)
library(mvtnorm)
library(flowEMMi)
```

    ## For detailed instructions please run browseVignettes('flowEMMi').
    ##   For an overview of available functions please run library(help='flowEMMi')

``` r
library(flowCore)
library(flowWorkspace)
```

    ## As part of improvements to flowWorkspace, some behavior of
    ## GatingSet objects has changed. For details, please read the section
    ## titled "The cytoframe and cytoset classes" in the package vignette:
    ## 
    ##   vignette("flowWorkspace-Introduction", "flowWorkspace")

``` r
library(ggcyto)
```

    ## Loading required package: ggplot2

    ## Loading required package: ncdfFlow

    ## Loading required package: BH

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ lubridate 1.9.2     ✔ tibble    3.2.1
    ## ✔ readr     2.1.4     ✔ tidyr     1.3.0

``` r
library(RColorBrewer)
library(knitr)
library(ellipse)
library(ggforce)
```

``` r
# create or *empty* the target directory, used to write this file's data: 
projthis::proj_create_dir_target(params$name, clean = TRUE)

# function to get path to target directory: path_target("sample.csv")
path_target <- projthis::proj_path_target(params$name)

# function to get path to previous data: path_source("00-import", "sample.csv")
path_source <- projthis::proj_path_source(params$name)
```

## Z-project

### Import preprocessed data

``` r
FDA_PI <- readRDS("~/Desktop/MSc_new_data/new_FDA_PI.rds")
gating_FDA_PI <- readRDS("~/Desktop/MSc_new_data/gating_FDA_PI.rds")
```

### flowEMMI gating on FDAPI

Based on the scatter plot in 01-data, we would like to set arrange,

- Setting range for x axis: 1000-60000

- Setting range for y axis: 1000-60000

``` r
set.seed(1)
location <- c("Inner_zone","Middle_zone","Outer_zone","Surrounding","Whole_colony")
gating_FDA_PI <- list()

for (i in 1:5){
  data_name <- paste0(location[i],"_FDA_PI.fcs")
  data <- FDA_PI[[i]]
  fdo <- mkFlowDataObject(data, xChannel="PMT.1", yChannel="PMT.3")
  gating <- flowEMMi( fdo=fdo, xMin=1000, xMax=60000, yMin=1000, yMax=60000
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
  gating_FDA_PI[[i]] <- gating$best
}
```

### Gating plots on FDA_PI

``` r
location <- c("Inner_zone","Middle_zone","Outer_zone","Whole_colony","Surrounding")
gating_FDA_PI_plot <- list()

for (i in 1:5){
  data_name <- paste0(location[i],"_FDA_PI.fcs")
  data <- FDA_PI[[i]]
  plots <- plotDensityAndEllipses(fcsData = data, ch1="PMT.1", ch2="PMT.3", alpha=0.9,
                            logScale = F, results = gating_FDA_PI[[i]],
                            title = data_name, plotRelevance = T,
                            ellipseDotSize = 0.5, axis_size=10, axisLabeling_size=10,
                            xlab = "Forward Scatter", ylab = "FDA_PI", font = "Arial")
  gating_FDA_PI_plot[[i]] <- plots$plot
}
```

    ## Warning in KernSmooth::bkde2D(x, bandwidth = bandwidth, gridsize = nbin, :
    ## Binning grid too coarse for current (small) bandwidth: consider increasing
    ## 'gridsize'

    ## Warning in KernSmooth::bkde2D(x, bandwidth = bandwidth, gridsize = nbin, :
    ## Binning grid too coarse for current (small) bandwidth: consider increasing
    ## 'gridsize'

![](02b-flowEMMI-on-FDA_PI-from-Z-project_files/figure-gfm/plot-FDAPI-1.png)<!-- -->![](02b-flowEMMI-on-FDA_PI-from-Z-project_files/figure-gfm/plot-FDAPI-2.png)<!-- -->![](02b-flowEMMI-on-FDA_PI-from-Z-project_files/figure-gfm/plot-FDAPI-3.png)<!-- -->

    ## Warning in KernSmooth::bkde2D(x, bandwidth = bandwidth, gridsize = nbin, :
    ## Binning grid too coarse for current (small) bandwidth: consider increasing
    ## 'gridsize'

![](02b-flowEMMI-on-FDA_PI-from-Z-project_files/figure-gfm/plot-FDAPI-4.png)<!-- -->![](02b-flowEMMI-on-FDA_PI-from-Z-project_files/figure-gfm/plot-FDAPI-5.png)<!-- -->

## Mahalanobis distance function with flowEMMi

Following is the integrated function.

Given a preprocessed dataset and the gating result from above, it can
generate the:

- Mahalanobis distance matrix & Cluster result

- Parameter table: Number of cells in each cluster, Area of the
  clustering ellipse, and the Coordinates of center

- Gating plots

``` r
flowEMMi_mahalanobis <- function(data,data_name,gating_data,alpha){
  
  mu <- gating_data@mu
  sigma <- gating_data@sigma
  
  names <- colnames(data)
  
  n_cells <- nrow(data)
  n_clusters <- length(sigma)
  
  #generate mahalanobis matrix
  maha_data <- matrix(NA,nrow=n_cells,ncol=n_clusters)
  
  for (i in 1:n_cells){
    for(j in 1:n_clusters){
      maha_data[i,j] <- mahalanobis(data[i,],mu[,j],sigma[[j]])
    }
  }
  
  maha <- maha_data[,2:n_clusters] %>% as.data.frame()
  
  #set 95% quantile as cutoff value
  threshold <- (-2*log(1-alpha))^2
  
  #determine cluster
  for (cell in 1:n_cells){
    rv <- maha[cell,1:n_clusters-1]
    if(all(rv>threshold)) { maha$Cluster[cell] <- NA}
    else { maha$Cluster[cell] <- which.min(rv)}
  }
  
  maha_data2 <- cbind(data,maha)
  
  test <- table(maha[,ncol(maha)]) %>% as.data.frame()
  coordinates <- sprintf("(%.2f,%.2f)",mu[1,2:ncol(mu)],mu[2,2:ncol(mu)])
  
  #Area of Ellipse
  eigen <- matrix(NA,nrow=length(sigma),ncol=2)
  for (i in 1:length(sigma)){
    eigen[i,] <- eigen(sigma[[i]])$values
  }
  eigen <- eigen[-1,]
  
  area <- matrix(NA,nrow=nrow(eigen),ncol=1)
  for (i in 1:nrow(eigen)){
    area[i,1] <- pi*sqrt(eigen[i,1]*eigen[i,2])
  }
  
  area <- area %>% as.data.frame()
  test <- cbind(test,area,coordinates)
  table <- test %>% 
    kable(caption = data_name,
          col.names = c("Cluster","Cells","Area","Coordinate"))
   
  #plot
  maha_data2$Cluster <- as.factor(maha_data2$Cluster)
  
  plot1 <- ggplot(maha_data2,aes(x=!!sym(names[1]),y=!!sym(names[2]),color=Cluster))+
    geom_point()+
    ggtitle(data_name)
  
  num_ellipse <- length(gating_data@sigma)
  
  for (j in 2:num_ellipse){
    mu <- gating_data@mu[,j]
    sigma <- gating_data@sigma[[j]]
    eli <- ellipse::ellipse(centre=mu,x=sigma,level=0.95,npoints=200) 
    eli <- as.data.frame(eli)
    colnames(eli)<- names
    plot1 <- plot1+geom_path(data = eli,
                  aes(x=!!sym(names[1]),y=!!sym(names[2])),color=j)
  }
  
  maha_result <- list(maha_matrix=maha_data2,table_info=table,plot=plot1)
  return(maha_result)
}
```

``` r
flowemmi_FDA_PI <- list()

for (i in 1:5){
  data <- FDA_PI[[i]]@exprs[,c(11,15)]
  data_name <- names(FDA_PI)[i]
  gating_data <- gating_FDA_PI[[i]]
  flowemmi_FDA_PI[[data_name]] <-flowEMMi_mahalanobis(data,data_name,gating_data,0.95)
  print(flowemmi_FDA_PI[[i]]$table_info)
  print(flowemmi_FDA_PI[[i]]$plot)
}
```

    ## 
    ## 
    ## Table: Inner_zone_FDA_PI.fcs
    ## 
    ## |Cluster | Cells|        Area|Coordinate          |
    ## |:-------|-----:|-----------:|:-------------------|
    ## |1       |  8644|    750346.6|(34342.30,54379.35) |
    ## |2       | 32205|   3020618.7|(34505.45,50486.47) |
    ## |3       |  2272|  10030909.7|(44862.67,31712.44) |
    ## |4       |  5060| 179762426.9|(21949.57,11431.48) |

![](02b-flowEMMI-on-FDA_PI-from-Z-project_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

    ## 
    ## 
    ## Table: Middle_zone_FDA_PI.fcs
    ## 
    ## |Cluster | Cells|     Area|Coordinate          |
    ## |:-------|-----:|--------:|:-------------------|
    ## |1       | 26756|  4694396|(35605.99,57083.86) |
    ## |2       |  2429| 11608940|(45918.96,26554.53) |
    ## |3       |  1220| 11631408|(33867.98,16704.77) |
    ## |4       |  1815| 24816958|(4941.35,8150.46)   |
    ## |5       | 17812| 19817643|(34366.51,33916.52) |

![](02b-flowEMMI-on-FDA_PI-from-Z-project_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

    ## 
    ## 
    ## Table: Outer_zone_FDA_PI.fcs
    ## 
    ## |Cluster | Cells|      Area|Coordinate          |
    ## |:-------|-----:|---------:|:-------------------|
    ## |1       | 44422|  18808340|(33281.33,31577.91) |
    ## |2       |  4329| 130033033|(15417.65,10506.33) |

![](02b-flowEMMI-on-FDA_PI-from-Z-project_files/figure-gfm/unnamed-chunk-2-3.png)<!-- -->

    ## 
    ## 
    ## Table: Surrounding_FDA_PI.fcs
    ## 
    ## |Cluster | Cells|     Area|Coordinate          |
    ## |:-------|-----:|--------:|:-------------------|
    ## |1       |  1440|  5540547|(46278.67,32603.67) |
    ## |2       |  3778|  4247893|(24430.99,4777.71)  |
    ## |3       |  3647|  4800805|(17510.04,2514.06)  |
    ## |4       |  8615| 12374104|(33950.58,13560.15) |

![](02b-flowEMMI-on-FDA_PI-from-Z-project_files/figure-gfm/unnamed-chunk-2-4.png)<!-- -->

    ## 
    ## 
    ## Table: Whole_colony_FDA_PI.fcs
    ## 
    ## |Cluster | Cells|      Area|Coordinate          |
    ## |:-------|-----:|---------:|:-------------------|
    ## |1       | 44621|  24276146|(34030.17,33927.28) |
    ## |2       |  1489|  11193828|(45646.15,30698.58) |
    ## |3       |  5608| 146142597|(17394.74,11857.46) |

![](02b-flowEMMI-on-FDA_PI-from-Z-project_files/figure-gfm/unnamed-chunk-2-5.png)<!-- -->

## Files written

These files have been written to the target directory,
`data/02b-flowEMMI on FDA_PI from Z-project`:

``` r
projthis::proj_dir_info(path_target())
```

    ## # A tibble: 0 × 4
    ## # ℹ 4 variables: path <fs::path>, type <fct>, size <fs::bytes>,
    ## #   modification_time <dttm>
