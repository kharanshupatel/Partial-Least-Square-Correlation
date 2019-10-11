---
title: 'Kharanshu Patel PLSC'
author: "Kharanshu Patel"
date: "10/24/2018"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

# devtools::install_github('HerveAbdi/PTCA4CATA', force = TRUE)
# devtools::install_github('HerveAbdi/data4PCCAR')
# install.packages("ggplot2")

library(ExPosition)
library(corrplot)
library(ggplot2)
library(InPosition)
library(dplyr)
library(PTCA4CATA)
library(data4PCCAR)
library(classInt)
library(Hmisc)
library(psych)
library(TInPosition)
library(RColorBrewer)
library(plotly)
library(heatmaply)
library(readxl)
library(TExPosition)

```
# Method: Partial Least Square Correlation (PLSC)

Partial Least Square (PLS) method compares the information in two data tables on the same set of observations. This method is also referred as projection to latent structures. To find the shared information between these two tables, the approach is equivalent to a correlation problem and the technique is then called partial least square correlation (PLSC). 
In PLSC, there are two sets of latent variables and these latent variables are required to have maximal covariance.
Latent variables are new variables obtained as linear combinations of the original variables.
These latent variables, which describe the observations, are required to “explain” the largest portion of the
covariance between the two tables. The original variables are described by their saliences.

# Dataset : oasis_cross-sectional

The Open Access Series of Imaging Studies (OASIS) is a project aimed at making MRI data sets of the brain freely available to the scientific community.
It is a cross-sectional MRI Data in Young, Middle Aged, Nondemented and Demented Older Adults.
This set consists of a cross-sectional collection of 216 subjects (rows) aged 33 to 96. The subjects include both men and women. 
It measures 216 subjects (rows) on 4 quanitiative variables (columns). 
For PLSC, the original dataset was divided into two separate tables based on the variables. Both tables consist of same rows.

```{r data_set}

# First: A clean start.
rm(list = ls())
graphics.off()

# read the file
oasis_cross_sectional <- read_excel("oasis_cross_sectional.xlsx")

# source the function file:
file4PLSCfunctions <- 'inferences4PLSC.R'
source(file4PLSCfunctions)

head(oasis_cross_sectional)
colnames(oasis_cross_sectional)[which(names(oasis_cross_sectional) == "M/F")] <- "Gender"

oasis_cross_sectional$Gender[oasis_cross_sectional$Gender =="M"] <- "1"
oasis_cross_sectional$Gender[oasis_cross_sectional$Gender =="F"] <- "2"

oasis_cross_sectional$Gender <- as.numeric(as.character(oasis_cross_sectional$Gender))

# Recode for PLSC
```{r PLSC recode}
## make 2 groups
## oasis.1 = MMSE, Age, Educ
## oasis.2 = CDR, eTiv, nWBV
## Design by Gender

oasis.1 <- oasis_cross_sectional[c(2,3,4)]
oasis.2 <- oasis_cross_sectional[c(6,8,9)]
oasis.des <- oasis_cross_sectional[,7]

## make a small scale example out of the data
## Here we take the first 150 observations
nN <- 216
o1 <- oasis.1[1:nN,]
o2 <- oasis.2[1:nN,]
od <- oasis.des[1:nN,]

```

# PLSC - DESIGN
```{r PLSC - DESIGN}

# Design colors according to CDR
DESIGN <- list()

DESIGN$CDR$vec <- od
DESIGN$CDR$mat <- makeNominalData(as.matrix(DESIGN$CDR$vec))
DESIGN$CDR$color_groups <- c("#e197ef", "#acd12b", "#687cef", "#750000")
DESIGN$CDR$color_obs <- as.matrix(DESIGN$CDR$vec)

DESIGN$CDR$color_obs[which(DESIGN$CDR$vec == "0")] <- 
  DESIGN$CDR$color_groups[1]
DESIGN$CDR$color_obs[which(DESIGN$CDR$vec == "0.5")] <- 
  DESIGN$CDR$color_groups[2]
DESIGN$CDR$color_obs[which(DESIGN$CDR$vec == "1")] <- 
  DESIGN$CDR$color_groups[3]
DESIGN$CDR$color_obs[which(DESIGN$CDR$vec == "2")] <- 
  DESIGN$CDR$color_groups[4]

# Colors for tables

col4fi <- c("#8b0000","#fd00e9","#fefc0e")
col4fj <- c("#006400","#00008b","#01fcff")

```

# PLSC
```{r PLSC}
# Run a PLSC with the epPLS() function
resPLSC <- tepPLS(o1,o2, DESIGN = od$CDR, graphs = FALSE)

eigs <- resPLSC$TExPosition.Data$eigs
tau <- resPLSC$TExPosition.Data$t
```

## PLSC - Permutation Test
```{r PLSC - Permutation Test}
resPerm4PLSC <- perm4PLSC(o1,  
                          o2, 
                          nIter = 1000
                          )

```

## PLSC - Bootstrap Ratios
```{r PLSC - Bootstrap Ratios}

resBoot4PLSC <- Boot4PLSC(o1, 
                          o2,
                          nIter = 1000,
                          Fi = resPLSC$TExPosition.Data$fi,
                          Fj = resPLSC$TExPosition.Data$fj,
                          nf2keep = 3,
                          eig = TRUE,
                          critical.value = 2,
                          alphaLevel = .05)

```

## PLSC - Correlation Plot
```{r PLSC - Heat Map Correlation}

corrplot::corrplot(cor(o1,o2),method = "circle",
         type = "full",
         order = "hclust",
         addrect = 2,
         col = c("white","black"), 
         bg = "gold2", 
         title = "PLSC Correlation Plot",
         mar=c(0,0,2,0), 
         tl.col = "#0c7315",
         tl.srt = 45,
         tl.cex = 1.3,
         na.label.col = "#fe7315",
         cl.cex = 0.8,
         cl.align.text = "l",)

```

# PLSC - Scree Plot + Interference
```{r PLSC - Scree Plot with Inference}

# PLSC - Scree Plot for inertia + Inference
PlotScree(ev = eigs, 
          p.ev =  resPerm4PLSC$pEigenvalues,
          title = 'PLSC Inertia Scree Plot with Inference',
          plotKaiser = TRUE
)

# PLSC - Scree Plot for singular values + Inference
PlotScree(ev = resPLSC$TExPosition.Data$eigs^(1/2),
          title = 'PLSC Singular Value Scree Plot with Inference',
          plotKaiser = TRUE, 
          color4Kaiser = ggplot2::alpha('darkorchid4', .5),
          lwd4Kaiser  = 2)

```

## PLSC - Salience (fi) for columns of table 1

```{r PLSC - Salience (fi), echo=TRUE}
fi <- resPLSC$TExPosition.Data$fi

sal.fi <- createFactorMap(fi,
                          axis1 = 1,
                          axis2 = 2,
                          display.points = TRUE,
                          display.names = TRUE,
                          display.labels = TRUE,
                          col.points = col4fi,
                          col.labels = col4fi,
                          alpha.points = 1,
                          pch = 20,
                          cex = 6,
                          text.cex = 5,
                          alpha.labels = 1,
                          col.axes = "black",
                          alpha.axes = 1,
                          col.background = "white",
                          width.axes = 1,
                          title = "PLSC Salience for Columns of Table 1")
# Labels for Inertia
label4Map.oasis <- createxyLabels.gen(1,2,
                                      lambda = eigs,
                                      tau = resPLSC$TExPosition.Data$t)

map.salFi <- sal.fi$zeMap + label4Map.oasis +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

print(map.salFi)

```

## PLSC - Salience (fj) for columns of table 2
```{r PLSC - Salience (fj), echo=TRUE}
fj <- resPLSC$TExPosition.Data$fj

sal.fj <- createFactorMap(fj,
                          display.points = TRUE,
                          display.names = TRUE,
                          display.labels = TRUE,
                          col.points = col4fj,
                          col.labels = col4fj,
                          alpha.points = 1,
                          pch = 20,
                          cex = 6,
                          text.cex = 5,
                          alpha.labels = 1,
                          col.axes = "black",
                          alpha.axes = 1,
                          col.background = "white",
                          width.axes = 1,
                          title = "PLSC Salience for Columns of Table 2")

# Labels for Inertia
label4Map.oasis <- createxyLabels.gen(1,2,
                                      lambda = eigs,
                                      tau = resPLSC$TExPosition.Data$t)

map.salFj <- sal.fj$zeMap + label4Map.oasis +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

dev.new()
print(map.salFj)

```

## PLSC - Latent variables for the rows from the first component
```{r PLSC - Latent}
Lx <- resPLSC$TExPosition.Data$lx
Ly <- resPLSC$TExPosition.Data$ly

laDim = 1
lv1.xy <- cbind(Lx[,laDim, drop = FALSE],
                Ly[,laDim, drop = FALSE])
colnames(lv1.xy) <- c(paste0('LX',laDim),paste0('LY',laDim))

# Latent Variable dimension 1
lat.1 <- createFactorMap(lv1.xy,
                         display.points = TRUE,
                         display.names = FALSE,
                         display.labels = FALSE,
                         col.points = DESIGN$CDR$color_obs,
                         col.labels = DESIGN$CDR$color_obs,
                         alpha.points = 0.5,
                         pch = 20,
                         cex = 2.8,
                         text.cex = 5,
                         alpha.labels = 1,
                         col.axes = "black",
                         alpha.axes = 0.6,
                         col.background = "white",
                         width.axes = 1,
                         title = "PLSC First Pair of Latent Variables : Lx1 + Ly1")

LV1.map <- lat.1$zeMap + xlab(paste0("X Latent Variable ", laDim)) +
  ylab(paste0("Y Latent Variable ", laDim))
print(LV1.map)

O.Means.1 <- PTCA4CATA::getMeans(lv1.xy, od$CDR)
rownames(O.Means.1) <- c("CDR-0", "CDR-0.5", "CDR-1", "CDR-2")

MapGroup.1 <- PTCA4CATA::createFactorMap(O.Means.1,
                                         display.labels = TRUE,
                                         col.points = DESIGN$CDR$color_groups,
                                         alpha.points = 1,
                                         cex = 6,  
                                         col.labels = "black",
                                         text.cex = 4,
                                         alpha.labels = 1,
                                         alpha.axes = 0.6,
                                         dev.new = FALSE,
                                         new.plot = FALSE)

latMeanPlot1 <- lat.1$zeMap + MapGroup.1$zeMap_dots + MapGroup.1$zeMap_text +
  ggtitle("PLSC Latent Variables : Lx1 + Ly1 with Means") + xlab(paste0("X Latent Variable ", laDim)) +
  ylab(paste0("Y Latent Variable ", laDim))

print(latMeanPlot1)
plsc.e1.l1 <- recordPlot()

# Latent Variable dimension 2
laDim = 2
lv2.xy <- cbind(Lx[,laDim, drop = FALSE],
                Ly[,laDim, drop = FALSE])
colnames(lv2.xy) <- c(paste0('LX',laDim),paste0('LY',laDim))

lat.2 <- createFactorMap(lv2.xy,
                         display.points = TRUE,
                         display.names = FALSE,
                         display.labels = FALSE,
                         col.points = DESIGN$CDR$color_obs,
                         col.labels = DESIGN$CDR$color_obs,
                         alpha.points = 0.5,
                         pch = 20,
                         cex = 2.8,
                         text.cex = 5,
                         alpha.labels = 1,
                         col.axes = "black",
                         alpha.axes = 0.6,
                         col.background = "white",
                         width.axes = 1,
                         title = "PLSC Second Pair of Latent Variables : Lx2 + Ly2")

LV2.map <- lat.2$zeMap + xlab(paste0("X Latent Variable ", laDim)) +
  ylab(paste0("Y Latent Variable ", laDim))
print(LV2.map)

O.Means.2 <- PTCA4CATA::getMeans(lv2.xy, od$CDR)
rownames(O.Means.1) <- c("CDR-0", "CDR-0.5", "CDR-1", "CDR-2")

MapGroup.2 <- PTCA4CATA::createFactorMap(O.Means.2,
                                         display.labels = TRUE,
                                         col.points = DESIGN$CDR$color_groups,
                                         alpha.points = 1,
                                         cex = 6,  
                                         col.labels = "black",
                                         text.cex = 4,
                                         alpha.labels = 1,
                                         alpha.axes = 0.6,
                                         dev.new = FALSE,
                                         new.plot = FALSE)

latMeanPlot2 <- lat.2$zeMap + MapGroup.2$zeMap_dots + MapGroup.2$zeMap_text +
  ggtitle("PLSC Latent Variables : Lx2 + Ly2 with Means") + xlab(paste0("X Latent Variable ", laDim)) +
  ylab(paste0("Y Latent Variable ", laDim))

print(latMeanPlot2)

# Bootstrap for CI:
BootCube.Gr.2 <- PTCA4CATA::Boot4Mean(lv2.xy, 
                                    design = oasis_cross_sectional$CDR,
                                    niter = 1000,
                                    suppressProgressBar = TRUE)

# Create Confidence Interval Plots
dimnames(BootCube.Gr.2$BootCube)[[2]] <- c("LX2","LY2")
GraphElli.2 <- PTCA4CATA::MakeCIEllipses(BootCube.Gr.2$BootCube[,1:2,],
                                         names.of.factors = c("LX2","LY2"),
                                         col = t(DESIGN$CDR$color_groups),
                                         p.level = .95
)
# create the I-map with Observations, means and confidence intervals
lv2.withCI <-  LV2.map + MapGroup.2$zeMap_text + 
  MapGroup.2$zeMap_dots +GraphElli.2

dev.new()
print(lv2.withCI)

```

# PLSC - Contributions
```{r PLSC - Contributions}

# Contributions for I
## Dimension 1
signed.ctrI <- resPLSC$TExPosition.Data$ci * sign(resPLSC$TExPosition.Data$fi)

b003.ctrI.s.1 <- PrettyBarPlot2(signed.ctrI[,1],
                                threshold = 1 / NROW(signed.ctrI),
                                font.size = 5,
                                color4bar = gplots::col2hex(col4fi), 
                                main = 'PLSC Contribution 1 for I',
                                ylab = 'Contributions',
                                ylim = c(1.2*min(signed.ctrI), 1.2*max(signed.ctrI)), horizontal = FALSE
)
print(b003.ctrI.s.1)
plsc.f1.c1.i <- recordPlot()

## Dimension 2
b003.ctrI.s.2 <- PrettyBarPlot2(signed.ctrI[,2],
                                threshold = 1 / NROW(signed.ctrI),
                                font.size = 5,
                                color4bar = gplots::col2hex(col4fi), 
                                main = 'PLSC Contribution 2 for I',
                                ylab = 'Contributions',
                                ylim = c(1.2*min(signed.ctrI), 1.2*max(signed.ctrI)), horizontal = FALSE
)
print(b003.ctrI.s.2)
plsc.f2.c2.i <- recordPlot()

# Contributions for J
## Dimension 1
signed.ctrJ <- resPLSC$TExPosition.Data$cj * sign(resPLSC$TExPosition.Data$fj)

b003.ctrJ.s.1 <- PrettyBarPlot2(signed.ctrJ[,1],
                                threshold = 1 / NROW(signed.ctrJ),
                                font.size = 5,
                                color4bar = gplots::col2hex(col4fj), 
                                main = 'PLSC Contribution 1 for J',
                                ylab = 'Contributions',
                                ylim = c(1.2*min(signed.ctrJ), 1.2*max(signed.ctrJ)), horizontal = FALSE
)
print(b003.ctrJ.s.1)
plsc.g1.c1.j <- recordPlot()

## Dimension 2
b003.ctrJ.s.2 <- PrettyBarPlot2(signed.ctrJ[,2],
                                threshold = 1 / NROW(signed.ctrJ),
                                font.size = 5,
                                color4bar = gplots::col2hex(col4fj), 
                                main = 'PLSC Contribution 2 for J',
                                ylab = 'Contributions',
                                ylim = c(1.2*min(signed.ctrJ), 1.2*max(signed.ctrJ)), horizontal = FALSE
)
print(b003.ctrJ.s.2)

```

# PLSC - Bootstrap
```{r PLSC - Bootstrap Plots}

BR.i <- resBoot4PLSC$bootRatios.i
BR.j <- resBoot4PLSC$bootRatios.j

BR.sig.i <- resBoot4PLSC$bootRatiosSignificant.i
BR.sig.j <- resBoot4PLSC$bootRatiosSignificant.j

# Bootstrap Ratios for I
## Dimension 1
laDim = 1
BR.i.map1 <- PrettyBarPlot2(BR.i[,laDim],
                            threshold = 2,
                            font.size = 5,
                            color4bar = gplots::col2hex(col4fi),
                            main = paste0(
                              'PLSC on the Oasis Set: Bootstrap ratio(I) Dimension ',laDim),
                            ylab = 'Bootstrap ratios'
)
print(BR.i.map1)
plsc.h1.br1.i <- recordPlot()

## Dimension 2
laDim = 2
BR.i.map2 <- PrettyBarPlot2(BR.i[,laDim],
                            threshold = 2,
                            font.size = 5,
                            color4bar = gplots::col2hex(col4fi),
                            main = paste0(
                              'PLSC on the Oasis Set: Bootstrap ratio(I) Dimension ',laDim),
                            ylab = 'Bootstrap ratios'
)
print(BR.i.map2)
plsc.h2.br2.i <- recordPlot()

# Bootstrap Ratios for J
## Dimension 1
laDim = 1
BR.j.map1 <- PrettyBarPlot2(BR.j[,laDim],
                            threshold = 2,
                            font.size = 5,
                            color4bar = gplots::col2hex(col4fj),
                            main = paste0(
                              'PLSC on the Oasis Set: Bootstrap ratio(J) Dimension ',laDim),
                            ylab = 'Bootstrap ratios'
)
print(BR.j.map1)
plsc.h3.br1.j <- recordPlot()

## Dimension 2
laDim = 2
BR.j.map2 <- PrettyBarPlot2(BR.j[,laDim],
                            threshold = 2,
                            font.size = 5,
                            color4bar = gplots::col2hex(col4fj),
                            main = paste0(
                              'PLSC on the Oasis Set: Bootstrap ratio(J) Dimension ',laDim),
                            ylab = 'Bootstrap ratios'
)
print(BR.j.map2)

```

# PLSC - Permutation Test for Eigen-values
```{r PLSC - Permutation Test for Eigen-values}

zeDim = 1
pH1 <- prettyHist(
  distribution = resPerm4PLSC$permEigenvalues[,zeDim], 
  observed = resPerm4PLSC$fixedEigenvalues[zeDim], 
  xlim = c(-0.5, 2), 
  breaks = 20,
  border = "white", 
  main = paste0("PLSC Permutation Test for Eigenvalue ",zeDim),
  xlab = paste0("Eigenvalue ",zeDim), 
  ylab = "", 
  counts = FALSE, 
  cutoffs = c( 0.975))
plsc.i1.perm1 <- recordPlot()

zeDim = 2
pH2 <- prettyHist(
  distribution = resPerm4PLSC$permEigenvalues[,zeDim], 
  observed = resPerm4PLSC$fixedEigenvalues[zeDim], 
  xlim = c(-0.2, 0.5),
  breaks = 20,
  border = "white", 
  main = paste0("PLSC Permutation Test for Eigenvalue ",zeDim),
  xlab = paste0("Eigenvalue ",zeDim), 
  ylab = "", 
  counts = FALSE, 
  cutoffs = c( 0.975))
plsc.i2.perm2 <- recordPlot()

zeDim = 3
pH3 <- prettyHist(
  distribution = resPerm4PLSC$permEigenvalues[,zeDim], 
  observed = resPerm4PLSC$fixedEigenvalues[zeDim], 
  xlim = c(-0.1, 0.1), 
  breaks = 20,
  border = "white", 
  main = paste0("PLSC Permutation Test for Eigenvalue ",zeDim),
  xlab = paste0("Eigenvalue ",zeDim), 
  ylab = "", 
  counts = FALSE, 
  cutoffs = c( 0.975))

```
# PLSC - Conclusion

Component 1
Salience:             Age   vs   nWBV
Latent Variable: 	 High CDR   vs   Low   CDR
Interpret: Higher chances of dementia with increase in age, also correspond to lesser whole brain volume.

Component 2
Salience:             Gender and eTIV
Latent Variable: 	 High CDR
Interpret: Females have higher chances of dementia and lower intracranial volume compared to Males.
