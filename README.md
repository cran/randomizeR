<!-- README.md is generated from README.Rmd. Please edit that file -->



# randomizeR

The R package **randomizeR** is an easy to use tool for the randomization of
clinical trials. Its aim is to provide sound scientific evidence for the choice
of a randomization procedure, on the basis of well-founded criteria. It includes
functions for generating randomization sequences from various randomization
procedures as well as various *issues* for the assessment of the randomization
procedures.

## Installation

You can install randomizeR from [CRAN](https://cran.r-project.org/package=randomizeR):

```r
install.packages('randomizeR', dependencies = TRUE)
```

## Motivation

Until now, the choice of a randomization procedure did not follow strict 
scientific principles. If the choice of a randomization procedure does not 
account for possible problems during the conduct of a clinical trial, the 
estimation of the treatment effect can be biased. randomizeR models the possible
problems that occur in the analysis of a clinical trial, and provides functions
to assess randomization procedures according to the issues.

## Usage

```r
library(randomizeR)
?randomizeR
vignette("comparison-example")

cg <- corGuess("CS")
rar <- genSeq(rarPar(6),10)
assess(rar, cg)
```

## License
The package is free and open source software, licensed under GPL 3.
