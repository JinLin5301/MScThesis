01-data
================
Compiled at 2023-08-13 20:54:45 UTC

``` r
here::i_am(paste0(params$name, ".Rmd"), uuid = "b3cd4969-aab1-4ff2-a553-9df4269240e4")
```

The purpose of this document is to introduce three flow cytometry data
that we’ll use for clustering comparison. They are:

-   The Z-project data

-   The microbial community data

-   The seaflow data

<!-- -->

    ## As part of improvements to flowWorkspace, some behavior of
    ## GatingSet objects has changed. For details, please read the section
    ## titled "The cytoframe and cytoset classes" in the package vignette:
    ## 
    ##   vignette("flowWorkspace-Introduction", "flowWorkspace")

    ## Loading required package: ggplot2

    ## Loading required package: ncdfFlow

    ## Loading required package: RcppArmadillo

    ## Loading required package: BH

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ lubridate 1.9.2     ✔ tibble    3.2.1
    ## ✔ readr     2.1.4     ✔ tidyr     1.3.0

``` r
# create or *empty* the target directory, used to write this file's data: 
projthis::proj_create_dir_target(params$name, clean = TRUE)

# function to get path to target directory: path_target("sample.csv")
path_target <- projthis::proj_path_target(params$name)

# function to get path to previous data: path_source("00-import", "sample.csv")
path_source <- projthis::proj_path_source(params$name)
```

## Z-project

Z-project is focusing on a colony of Bacillus subtilis, lab exams the
sub-population structure at five different locations: inner zone, middle
zone, outer zone, surrounding and whole colony.

``` r
# import data
fcs.direction <- "~/Desktop/Z-Project Bacillus_pretest/fcs/"
fcs.file.name <- list.files(fcs.direction,pattern = "\\.fcs$",full.names = T)
z_fs <- list()

for (i in fcs.file.name){
  data <- read.FCS(i,alter.names = T,transformation = F)
  new_name <- substr(i,start = 55,stop=nchar(i))
  z_fs[[new_name]] <- data
}

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

DAPI,FDA,PI are common fluorescent dyes which are excited with 488nm and
measured with 460/50nm, 530/40nm, 616/23nm respectively.

Based on above table, three variables for these three fluorescent dyes
are: DAPI-PMT.9, FDA-PMT.3, PI-PMT.5.

Now, let’s have a look of cytograms from Z-project.

``` r
plot_DAPI_list <- list()
plot_FDA_PI_list <- list()

for (i in 1:5){
  data_name <- paste0(location[i],"_DAPI.fcs")
  data <- new_DAPI_fs[[data_name]]
  p <- ggcyto(data, aes(x = "PMT.1", y = "PMT.9")) + geom_hex(bins = 128)
  plot_DAPI_list[[i]] <- p
}

for (i in 1:5){
  data_name <- paste0(location[i],"_FDA_PI.fcs")
  data <- new_FDA_PI_fs[[data_name]]
  p <- ggcyto(data, aes(x = "PMT.1", y = "PMT.3")) + geom_hex(bins = 128)
  plot_FDA_PI_list[[i]] <- p
}

grid.arrange(grobs = plot_DAPI_list, ncol = 3)
grid.arrange(grobs = plot_FDA_PI_list, ncol = 3)
```

This is another visulization method offered by flowEMMI v2 packages.

``` r
# replace the fcsData and ch1/2
p <- plotDensityAndEllipses(fcsData = data, ch1="PMT.1", ch2="PMT.9",
                          logScale = F, title = data_name,
                          axis_size=10, axisLabeling_size=10,
                          xlab = "Forward Scatter", ylab = "DAPI Fluorescence", 
                          font = "Arial")
```

## Microbial community

For this section, we have two datasets. The data from Liu et. al doesn’t
have mass transfer but it contains a abiotic parameter measurement.
Another is from Li et. al and by contrast, it has looped mass transfter
but no environment covariable records.

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
reactor <- c("C1","C2","D1","D2","D3")
data_list <- list()

# get file path for each reactor
for (i in reactor){
  path <- "~/Desktop/FlowRepository_Liu/"
  path <- paste0(path,i,"/")
  file_name <- list.files(path,pattern = ".fcs$",full.names = T)
  for (j in file_name){
    data <- read.FCS(j,alter.names=T,transformation=F)
    new_name <- substr(j,start=44,stop=58) # use nchar
    data_list[[new_name]] <- data
  }
}

# remove the debris
new_data <- list()

n_data <- length(data_list)

for (i in 1:n_data) {
  data_name <- names(data_list)[i]
  data <- data_list[[i]]
  sub <- data@exprs %>% as.data.frame() %>% dplyr::filter(PMT.1>200) %>% dplyr::filter(PMT.9>1000)
  sub <- as.matrix(sub)
  data@exprs <- sub
  new_data[[data_name]]<-data
}

# optical data
abio_data <- read_excel("~/Desktop/osfstorage/Data sheet for abiotic parameters.xlsx")
```

### Li et. al

Five local microbiomes were grown for over 114 generations and connected
by a loop. The flow cytometry raw data has been deposited in
FlowRepository (<https://flowrepository.org/>) with accession number:
FR-FCM-Z3MU

``` r
folder_path <- "~/Desktop/FlowRepo/"
file_name <- list.files(folder_path,pattern = "\\.fcs$",full.names = T)
MoFloA_Feb_data <- list()
MoFloA_Mar_data <- list()
MoFloA_Nov_data <- list()
MoFloA_Dec_data <- list()

# Group files by the month

for (file_name in file_name){
  month <- substr(file_name,start = 46,stop = 48) # use nchar to locate
  new_name <- substr(file_name,start = 43, stop = 58)
  data <- read.FCS(file_name,alter.names = T,transformation=F)
  if (month=="Feb") { MoFloA_Feb_data[[new_name]] <- data}
  else if (month=="Mar") {MoFloA_Mar_data[[new_name]] <- data}
  else if (month=="Nov") {MoFloA_Nov_data[[new_name]] <- data}
  else if (month=="Dec") {MoFloA_Dec_data[[new_name]] <- data}
}
```

## Seaflow

``` r
seaflow <- readRDS(file="~/Desktop/paper-data/MGL1704-hourly-paper.RDS")
```

## Files written

These files have been written to the target directory, data/01-data:

``` r
projthis::proj_dir_info(path_target())
```

    ## # A tibble: 0 × 4
    ## # ℹ 4 variables: path <fs::path>, type <fct>, size <fs::bytes>,
    ## #   modification_time <dttm>
