01-data
================
Compiled at 2023-09-20 11:12:21 UTC

The purpose of this document is to introduce three flow cytometry data
that we’ll use for clustering comparison. They are:

- The Z-project data

- The microbial community data

- The seaflow data

``` r
library("conflicted")
library(readxl)
library(purrr)
library(dplyr)
library(flowCore)
#library(flowmix)
library(flowWorkspace)
library(ggcyto)
library(tidyverse)
library(gridExtra)
library(knitr)
library(ggplot2)
library(RColorBrewer)
```

## Z-project

Z-project is focusing on a colony of Bacillus subtilis, lab exams the
sub-population structure at five different locations: inner zone, middle
zone, outer zone, surrounding and whole colony.

``` r
# import data
fcs.direction <- "~/MSc_raw_data/Z-Project Bacillus_pretest/fcs/"
fcs.file.name <- list.files(fcs.direction,pattern = "\\.fcs$",full.names = T)
z_fs <- list()

for (i in fcs.file.name){
  data <- read.FCS(i,alter.names = T,transformation = F)
  new_name <- substr(i,start = 72,stop=nchar(i))
  z_fs[[new_name]] <- data
}
```

``` r
markernames(z_fs[[2]])
```

    ##          PMT.1        X.PMT.1          PMT.2        X.PMT.2          PMT.3 
    ##    "FSC [488]"    "FSC [488]"    "SSC [488]"    "SSC [488]" "530/40 [488]" 
    ##        X.PMT.3          PMT.4        X.PMT.4          PMT.5        X.PMT.5 
    ## "530/40 [488]" "580/30 [488]" "580/30 [488]" "616/23 [488]" "616/23 [488]" 
    ##          PMT.6        X.PMT.6          PMT.7        X.PMT.7          PMT.8 
    ## "780/60 [488]" "780/60 [488]" "670/30 [640]" "670/30 [640]" "730/45 [640]" 
    ##        X.PMT.8          PMT.9        X.PMT.9         PMT.10       X.PMT.10 
    ## "730/45 [640]" "460/50 [355]" "460/50 [355]"  "650LP [355]"  "650LP [355]" 
    ##         PMT.11       X.PMT.11         PMT.12       X.PMT.12         PMT.13 
    ## "460/50 [405]" "460/50 [405]" "520/35 [405]" "520/35 [405]" "605/40 [405]" 
    ##       X.PMT.13         PMT.14       X.PMT.14     PMT.1.Area    PMT.1.Width 
    ## "605/40 [405]"  "650LP [405]"  "650LP [405]"    "FSC [488]"    "FSC [488]" 
    ##     PMT.2.Area    PMT.2.Width 
    ##    "SSC [488]"    "SSC [488]"

DAPI, FDA, PI are common fluorescent dyes which are excited with 488nm
and measured with 460/50nm, 530/40nm, 616/23nm respectively.

Based on above table, the relationship between fluorescent dyes and
variables are: DAPI-PMT.9, FDA-PMT.3, PI-PMT.5.

``` r
location <- c("Inner_zone","Middle_zone","Outer_zone","Surrounding","Whole_colony")

# remove the debris
new_DAPI_fs <- list()
new_FDA_PI_fs <- list()

for (i in 1:5) {
  data_name <- paste0(location[i],"_DAPI.fcs")
  data <- z_fs[[data_name]]
  sub <- data@exprs %>% as.data.frame() %>% dplyr::filter(PMT.1>200) %>% dplyr::filter(PMT.9>200)
  sub <- as.matrix(sub)
  data@exprs <- sub
  new_DAPI_fs[[data_name]]<-data
}

for (i in 1:5) {
  data_name <- paste0(location[i],"_FDA_PI.fcs")
  data <- z_fs[[data_name]]
  sub <- data@exprs %>% as.data.frame() %>% dplyr::filter(PMT.1>200) %>% dplyr::filter(PMT.3>200) %>% dplyr::filter(PMT.5>200)
  sub <- as.matrix(sub)
  data@exprs <- sub
  new_FDA_PI_fs[[data_name]]<-data
}
```

``` r
write.FCS(new_DAPI_fs,"~/MScThesis/workflow/data/01-data/new_DAPI.fcs")
write.FCS(new_FDA_PI_fs,"~/MScThesis/workflow/data/01-data/new_FDA_PI.fcs")
```

Now, let’s have a look of cytograms from Z-project.

``` r
for (i in 1:5){
  data_name <- paste0(location[i],"_DAPI.fcs")
  data <- new_DAPI_fs[[data_name]]
  p <- ggcyto(data, aes(x = "PMT.1", y = "PMT.9")) + geom_hex(bins = 128)
  print(p)
}
```

![](01-data_files/figure-gfm/plot-z-1.png)<!-- -->![](01-data_files/figure-gfm/plot-z-2.png)<!-- -->![](01-data_files/figure-gfm/plot-z-3.png)<!-- -->![](01-data_files/figure-gfm/plot-z-4.png)<!-- -->![](01-data_files/figure-gfm/plot-z-5.png)<!-- -->

``` r
for (i in 1:5){
  data_name <- paste0(location[i],"_FDA_PI.fcs")
  data <- new_FDA_PI_fs[[data_name]]
  p <- ggcyto(data, aes(x = "PMT.1", y = "PMT.3")) + geom_hex(bins = 128)
  print(p)
}
```

![](01-data_files/figure-gfm/plot-z-6.png)<!-- -->![](01-data_files/figure-gfm/plot-z-7.png)<!-- -->![](01-data_files/figure-gfm/plot-z-8.png)<!-- -->![](01-data_files/figure-gfm/plot-z-9.png)<!-- -->![](01-data_files/figure-gfm/plot-z-10.png)<!-- -->

This is another visulization method offered by flowEMMI v2 packages.

``` r
# replace the fcsData and ch1/2
p <- plotDensityAndEllipses(fcsData = data, ch1="PMT.1", ch2="PMT.x",
                          logScale = F, title = data_name,
                          axis_size=10, axisLabeling_size=10,
                          xlab = "Forward Scatter", ylab = "Fluorescence", 
                          font = "Arial")
```

## Microbial community

For this section, we have two datasets. The data from Liu et. al doesn’t
have mass transfer but it contains a abiotic parameter measurement.
Another is from Li et. al, by contrast, it has looped mass transfter but
no environment covariable records.

### Liu et. al

This flow cytometry data consists of 5 reactors data.

C1 and C2 are two control reactors in order to test for neutral
behaviour under undisturbed communities.

D1, D2 and D3 were disturbed by a repeated and soft temperature stressor
from 30 to 40 degrees.

Liu et. al sampled from all five reactors and analyzed them by flow
cytometry. All raw data can be accessed by flowRepository
(<https://flowrepository.org/>) under accession number: FR-FCM-ZYWX

``` r
# reactor <- c("C1","C2","D1","D2","D3")
reactor <- c("C1")
micro_data_list <- list()

# get file path for each reactor
for (i in reactor){
  path <- "~/MSc_raw_data/FlowRepository_Liu/"
  path <- paste0(path,i,"/")
  file_name <- list.files(path,pattern = ".fcs$",full.names = T)
  for (j in file_name){
    data <- read.FCS(j,alter.names=T,transformation=F)
    new_name <- substr(j,start=63,stop=77) # use nchar 
    micro_data_list[[new_name]] <- data
  }
}

# optical data
abio_data <- read_excel("~/MSc_raw_data/Data sheet for abiotic parameters.xlsx")
```

``` r
markernames(micro_data_list[[1]])
```

    ##          PMT.1        X.PMT.1          PMT.2        X.PMT.2          PMT.3 
    ##    "FSC [488]"    "FSC [488]"    "SSC [488]"    "SSC [488]" "530/40 [488]" 
    ##        X.PMT.3          PMT.4        X.PMT.4          PMT.5        X.PMT.5 
    ## "530/40 [488]" "580/30 [488]" "580/30 [488]" "616/23 [488]" "616/23 [488]" 
    ##          PMT.6        X.PMT.6          PMT.7        X.PMT.7          PMT.8 
    ## "780/60 [488]" "780/60 [488]" "670/30 [640]" "670/30 [640]" "730/45 [640]" 
    ##        X.PMT.8          PMT.9        X.PMT.9         PMT.10       X.PMT.10 
    ## "730/45 [640]" "460/50 [355]" "460/50 [355]"  "650LP [355]"  "650LP [355]" 
    ##         PMT.11       X.PMT.11         PMT.12       X.PMT.12         PMT.13 
    ## "460/50 [405]" "460/50 [405]" "520/35 [405]" "520/35 [405]" "605/40 [405]" 
    ##       X.PMT.13         PMT.14       X.PMT.14 
    ## "605/40 [405]"  "650LP [405]"  "650LP [405]"

In this case, we’re interested in the forwards scatter, indicated the
cell size, named as PMT.1, and fluorescence DAPI, named as PMT.9.

So we could remove the debris now.

``` r
# remove the debris
new_micro_data <- list()

n_data <- length(micro_data_list)

for (i in 1:n_data) {
  data_name <- names(micro_data_list)[i]
  data <- micro_data_list[[i]]
  sub <- data@exprs %>% as.data.frame() %>% dplyr::filter(PMT.1>200) %>% dplyr::filter(PMT.9>1000)
  sub <- as.matrix(sub)
  data@exprs <- sub
  new_micro_data[[data_name]]<-data
}
```

After removing the debris, we’d like to separate the whole dataset into
five subset, represented by the reactor names.

``` r
# separate by reactor
C1_data <- list()
C2_data <- list()
D1_data <- list()
D2_data <- list()
D3_data <- list()

for (i in names(micro_data_list)){
  data <- new_micro_data[[i]]
  idx <- substr(i,start = 10,stop=11)
  if(idx=="C1"){C1_data[[i]]<-data}
  else if(idx=="C2"){C2_data[[i]]<-data}
  else if(idx=="D1"){D1_data[[i]]<-data}
  else if(idx=="D2"){D2_data[[i]]<-data}
  else if(idx=="D3"){D3_data[[i]]<-data}
}
```

``` r
for (i in 1:65){
  data <- C1_data[[i]]
  p <- ggcyto(data, aes(x = "PMT.1", y = "PMT.9")) + geom_hex(bins = 128)
  print(p)
}
```

![](01-data_files/figure-gfm/plot-micro-1.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-2.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-3.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-4.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-5.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-6.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-7.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-8.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-9.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-10.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-11.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-12.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-13.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-14.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-15.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-16.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-17.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-18.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-19.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-20.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-21.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-22.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-23.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-24.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-25.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-26.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-27.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-28.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-29.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-30.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-31.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-32.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-33.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-34.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-35.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-36.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-37.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-38.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-39.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-40.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-41.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-42.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-43.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-44.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-45.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-46.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-47.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-48.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-49.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-50.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-51.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-52.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-53.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-54.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-55.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-56.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-57.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-58.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-59.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-60.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-61.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-62.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-63.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-64.png)<!-- -->![](01-data_files/figure-gfm/plot-micro-65.png)<!-- -->

The table reveals that this dataset is too costly for the downstream
analysis, so we’d like to bin the data first to reduce the complexity of
raw data.

Before doing the binning, we need to transfer the flowFrame to list
first.

``` r
set.seed(0)
micro_data_ylist <- list()
grid_micro_ylist <- list()
grid_micro_countslist <- list()

for (i in seq_along(new_micro_data)){
  micro_data_ylist[[i]] <- new_micro_data[[i]]@exprs[,c(11,13,27)] # PMT.1/2 for FSC/SSC, PMT.9 for fluorescence
  micro_data_ylist[[i]][,] <- as.double(micro_data_ylist[[i]][,])
  grid <- make_grid(micro_data_ylist[[i]],gridsize = 30)
  obj <- bin_many_cytograms(micro_data_ylist,grid,mc.cores = 4,verbose = F)
  grid_micro_ylist[[i]] <- obj$ybin_list
  grid_micro_countslist[[i]] <- obj$counts_list
}
```

``` r
write.FCS(new_micro_data,"~/MScThesis/workflow/data/01-data/new_micro_data.fcs")
```

### Li et. al

Five local microbiomes were grown for over 114 generations and connected
by a loop. The flow cytometry raw data has been deposited in
FlowRepository (<https://flowrepository.org/>) with accession number:
FR-FCM-Z3MU

``` r
folder_path <- "~/MSc_raw_data/FlowRepo_Li/"
file_name <- list.files(folder_path,pattern = "\\.fcs$",full.names = T)
MoFloA_Feb_data <- list()
MoFloA_Mar_data <- list()
MoFloA_Nov_data <- list()
MoFloA_Dec_data <- list()

# group files by the month
for (file_name in file_name){
  month <- substr(file_name,start = 66,stop = 68) # use nchar to locate 
  new_name <- substr(file_name,start = 53, stop = 77)
  data <- read.FCS(file_name,alter.names = T,transformation=F)
  if (month=="Feb") { MoFloA_Feb_data[[new_name]] <- data}
  #else if (month=="Mar") {MoFloA_Mar_data[[new_name]] <- data}
  #else if (month=="Nov") {MoFloA_Nov_data[[new_name]] <- data}
  #else if (month=="Dec") {MoFloA_Dec_data[[new_name]] <- data}
}
```

From the instrumental setup, we could tell that the DAPI fluorescenece
was measured at signal channel FL4.

``` r
markernames(MoFloA_Feb_data[[21]])
```

    ##       Pulse.Width            FS.Lin            SS.Lin          FL.1.Lin 
    ##     "Pulse Width"              "FS"              "SS"            "FL 1" 
    ##          FL.2.Lin          FL.3.Lin          FL.4.Lin          FL.5.Lin 
    ##            "FL 2"            "FL 3"            "FL 4"            "FL 5" 
    ##          FL.6.Lin   Sort.Classifier 
    ##            "FL 6" "Sort Classifier"

This is an example that the variable was wrongly recorded. Therefore,
before further analysis, we’d like to correct these mistakes.

``` r
# remove the debris
new_MoFloA_Feb <- list()
Q_data <- c()

for (i in 1:268) {
  data_name <- names(MoFloA_Feb_data)[i]
  data1 <- MoFloA_Feb_data[[i]]
  annodata <- data1@parameters@data
  if(any(grepl("Lin",annodata[[1]]))){
    Q_data <- c(Q_data,i)
  }
}

print(paste0("This record needs to be removed: ",Q_data))
```

    ## [1] "This record needs to be removed: 21" "This record needs to be removed: 31"
    ## [3] "This record needs to be removed: 61"

``` r
new_MoFloA_Feb_data <- MoFloA_Feb_data[-Q_data]

# removing debris
for (i in 1:length(new_MoFloA_Feb_data)){
  data_name <- names(new_MoFloA_Feb_data)[i]
  data <- new_MoFloA_Feb_data[[i]]
  sub <- data@exprs %>% as.data.frame() 
  sub <- sub %>% dplyr::filter(FL.4.Log>200) 
  sub <- as.matrix(sub)
  data@exprs <- sub
  new_MoFloA_Feb[[data_name]]<-data
}
```

``` r
for (i in 1:25){
  data <- new_MoFloA_Feb[[i]]
  p <- ggcyto(data, aes(x = `FS.Log`, y = `FL.4.Log`)) + geom_hex(bins = 128)
  print(p)
}
```

![](01-data_files/figure-gfm/plot-micro2-1.png)<!-- -->![](01-data_files/figure-gfm/plot-micro2-2.png)<!-- -->![](01-data_files/figure-gfm/plot-micro2-3.png)<!-- -->![](01-data_files/figure-gfm/plot-micro2-4.png)<!-- -->![](01-data_files/figure-gfm/plot-micro2-5.png)<!-- -->![](01-data_files/figure-gfm/plot-micro2-6.png)<!-- -->![](01-data_files/figure-gfm/plot-micro2-7.png)<!-- -->![](01-data_files/figure-gfm/plot-micro2-8.png)<!-- -->![](01-data_files/figure-gfm/plot-micro2-9.png)<!-- -->![](01-data_files/figure-gfm/plot-micro2-10.png)<!-- -->![](01-data_files/figure-gfm/plot-micro2-11.png)<!-- -->![](01-data_files/figure-gfm/plot-micro2-12.png)<!-- -->![](01-data_files/figure-gfm/plot-micro2-13.png)<!-- -->![](01-data_files/figure-gfm/plot-micro2-14.png)<!-- -->![](01-data_files/figure-gfm/plot-micro2-15.png)<!-- -->![](01-data_files/figure-gfm/plot-micro2-16.png)<!-- -->![](01-data_files/figure-gfm/plot-micro2-17.png)<!-- -->![](01-data_files/figure-gfm/plot-micro2-18.png)<!-- -->![](01-data_files/figure-gfm/plot-micro2-19.png)<!-- -->![](01-data_files/figure-gfm/plot-micro2-20.png)<!-- -->![](01-data_files/figure-gfm/plot-micro2-21.png)<!-- -->![](01-data_files/figure-gfm/plot-micro2-22.png)<!-- -->![](01-data_files/figure-gfm/plot-micro2-23.png)<!-- -->![](01-data_files/figure-gfm/plot-micro2-24.png)<!-- -->![](01-data_files/figure-gfm/plot-micro2-25.png)<!-- -->

``` r
write.FCS(new_MoFloA_Feb,"~/MScThesis/workflow/data/01-data/new_MoFloA_Feb.fcs")
```

## Seaflow

``` r
seaflow <- readRDS(file="~/MSc_raw_data/MGL1704-hourly-paper.RDS")
seaflow_ylist <- seaflow$ylist
seaflow_X <- seaflow$X

print(paste("This dataset contains",length(seaflow_ylist),"time points."))
```

    ## [1] "This dataset contains 296 time points."

``` r
print(paste("There are",ncol(seaflow_X),"environment variables."))
```

    ## [1] "There are 39 environment variables."

Since this dataset has already binned into a list, former visualization
methods for fcs files aren’t suitable anymore. So, we apply the basic
scatter plots here to illustrate the relationship between cell size and
the fluorescence.

``` r
for ( i in 1:20){
  data <- seaflow_ylist[[i]] %>% as.data.frame()
  data_name <- names(seaflow_ylist)[i]
  p <- ggplot(data,aes(x=diam_mid,y=chl_small))+geom_hex(bins=32)+
    scale_fill_gradient(low="blue",high="red")+
    labs(x="FSC",y="Red fluorescence",title=data_name)
  print(p)
}
```

![](01-data_files/figure-gfm/plot-sea-1.png)<!-- -->![](01-data_files/figure-gfm/plot-sea-2.png)<!-- -->![](01-data_files/figure-gfm/plot-sea-3.png)<!-- -->![](01-data_files/figure-gfm/plot-sea-4.png)<!-- -->![](01-data_files/figure-gfm/plot-sea-5.png)<!-- -->![](01-data_files/figure-gfm/plot-sea-6.png)<!-- -->![](01-data_files/figure-gfm/plot-sea-7.png)<!-- -->![](01-data_files/figure-gfm/plot-sea-8.png)<!-- -->![](01-data_files/figure-gfm/plot-sea-9.png)<!-- -->![](01-data_files/figure-gfm/plot-sea-10.png)<!-- -->![](01-data_files/figure-gfm/plot-sea-11.png)<!-- -->![](01-data_files/figure-gfm/plot-sea-12.png)<!-- -->![](01-data_files/figure-gfm/plot-sea-13.png)<!-- -->![](01-data_files/figure-gfm/plot-sea-14.png)<!-- -->![](01-data_files/figure-gfm/plot-sea-15.png)<!-- -->![](01-data_files/figure-gfm/plot-sea-16.png)<!-- -->![](01-data_files/figure-gfm/plot-sea-17.png)<!-- -->![](01-data_files/figure-gfm/plot-sea-18.png)<!-- -->![](01-data_files/figure-gfm/plot-sea-19.png)<!-- -->![](01-data_files/figure-gfm/plot-sea-20.png)<!-- -->

In order to apply this dataset on other flow cytometry methods, we’d
like to convert it into flowFrame object.

``` r
template <- z_fs[[1]]
yname <- names(seaflow_ylist)

sea_flowframe <- list()

for (i in 1:296){
  sea_flowframe <- c(sea_flowframe,template)
  names(sea_flowframe)[i]<- yname[i]
  sea_flowframe[[i]]@exprs <- seaflow_ylist[[i]]
  annodata <- template@parameters@data[1:3,]
  annodata[,1] <- c("diam_mid","chl_small","pe")
  annodata[,2] <- c("diameter","prochlorococcus","<NA>")
  annodata[,3] <- c(8,8,8)
  annodata[,5] <- c(8,8,8)
  sea_flowframe[[i]]@parameters@data <- annodata
}
```

``` r
for (i in 1:10){
  data_name <- names(sea_flowframe)[i]
  data <- sea_flowframe[[data_name]]
  p <- ggcyto(data, aes(x = "diam_mid", y = "chl_small")) + geom_hex(bins = 128)+ggtitle(data_name)
  print(p)
}
```

![](01-data_files/figure-gfm/plot-seaframe-1.png)<!-- -->![](01-data_files/figure-gfm/plot-seaframe-2.png)<!-- -->![](01-data_files/figure-gfm/plot-seaframe-3.png)<!-- -->![](01-data_files/figure-gfm/plot-seaframe-4.png)<!-- -->![](01-data_files/figure-gfm/plot-seaframe-5.png)<!-- -->![](01-data_files/figure-gfm/plot-seaframe-6.png)<!-- -->![](01-data_files/figure-gfm/plot-seaframe-7.png)<!-- -->![](01-data_files/figure-gfm/plot-seaframe-8.png)<!-- -->![](01-data_files/figure-gfm/plot-seaframe-9.png)<!-- -->![](01-data_files/figure-gfm/plot-seaframe-10.png)<!-- -->

``` r
write.FCS(sea_flowframe,"~/MScThesis/workflow/data/01-data/sea_flowframe.fcs")
```

## Files written

These files have been written to the target directory, data/01-data:

``` r
projthis::proj_dir_info(path_target())
```

    ## # A tibble: 0 × 4
    ## # ℹ 4 variables: path <fs::path>, type <fct>, size <fs::bytes>,
    ## #   modification_time <dttm>
