# ShrinkBayes
Updated 2-2-2012: ShrinkBayes can now deal with different designs (and different number of data points) per feature

R package for differential expression: sequencing counts (RNAseq), microarrays, and HTRNAi

This package corresponds to the papers:

- Van de Wiel MA, Neerincx M, Buffart TE, Sie D, Verheul HMW (2014). ShrinkBayes: a versatile R-package for analysis of count-based sequencing data in complex study designs. BMC Bioinformatics. 15(1):116.

- Van de Wiel MA,  De Menezes RX, Siebring E, Van Beusechem VW (2013). Analysis of small-sample clinical genomics studies using multi-parameter shrinkage: application to high-throughput RNA interference screening. BMC Med Genom 6 (Suppl 2), S1.

- Van de Wiel MA, Leday GGR, Pardo L, Rue H, Van der Vaart AW, Van Wieringen WN (2012). Bayesian analysis of RNA sequencing data by estimating multiple shrinkage priors. Biostatistics 14, 113-128.


Why use ShrinkBayes?
•It is demonstrated to be more reproducible than most other approaches, in particular for small samples
•It is build upon the INLA fundament. Hence, much faster than MCMC.
•Allows for (zero-inflated) Counts and Gaussian data [Hence can be applied to (mi)RNAseq, CAGE, mRNA/miRNA microaarray, HT RNAi, etc....]

•Very versatile in terms of designs. GLM-context, and allows for random effects.
•Provides Bayesian FDR and lfdr estimates.
•Accommodates a variety of priors, including mixture and nonparametric ones.
•Enables multi-parameter shrinkage

R Package
Note: if you have a choice to use either Windows or Unix/Linux, opt for the latter. ShrinkBayes runs more efficiently under Unix/Linux than under Windows. NOTE:  when running ShrinkBayes you may see *** WARNINGS ***  from INLA (e.g. on eigenvalues, or on convergence, or even something like 18500 Aborted...). They can currently not be surpressed, because they are produced by C-code. Please ignore them. 

Installation instructions
 For Windows users: PLEASE shutdown Windows Error Reporting. Windows XP: Windows key + Pause/Break, Advanced, Error reporting, Completely (no critical errors either). Windows 7 (and other) see: shutdown error reporting

ShrinkBayes depends on the following packages (see below for installation): 
INLA  (which requires packages sp and pixmap), snowfall, VGAM, mclust, logcondens, Iso, XML, rgl [All available from CRAN]

Steps:
 1. install.packages(c("sp","pixmap", "snowfall", "VGAM", "mclust", "logcondens", "Iso","XML","rgl"), repos="http://cran.r-project.org")

Unix/Linux: if you can't install "XML", "rgl", try
sudo apt-get build-dep r-cran-xml
sudo apt-get build-dep r-cran-rgl

2.source("http://www.math.ntnu.no/inla/givemeINLA.R") 
 [or if you installed INLA before 01/10/2012 you should upgrade by using inla.upgrade() ]

Then, install ShrinkBayes:
library("devtools")
install_github("markvdwiel/ShrinkBayes")

A pdf Vignette is available in the inst/doc folder
