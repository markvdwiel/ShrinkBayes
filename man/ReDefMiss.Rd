\name{ReDefMiss}
\alias{ReDefMiss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Redefines missing covariates.
}
\description{
Depending on the type of covariate (factor, numeric binary, continuous) and the interpretation of a missing, 
this function assigns a value to missings and possibly rescales the continuous variable or recodes 
a numeric binary one.
}
\usage{
ReDefMiss(para, NAzero = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{para}{
Character. Name of the variable.
}
  \item{NAzero}{
Boolean. If \code{TRUE} NAs are simply replaced by 0; other values are left unchanged.
}
}
\details{
If \code{NAzero=FALSE}: then if the variable is a numeric binary one (0,1 values), 0's are replaced by -1 and 
NAs are replaced by 0; if the variable is a factor, NAs are replaced by level "0" of the factor; 
if the variable is continuous (or ratio-scaled) then the median of the non-missing values is subtracted from those 
and NAs are replace by 0. 
}
\value{
Factor or numeric vector, representing the recoded variable
}

\author{
Mark A. van de Wiel
}



