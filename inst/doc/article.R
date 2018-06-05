## ----setup, include = FALSE, cache = FALSE--------------------------
library(knitr)
# set global chunk options
knitr::render_sweave() 
#opts_chunk$set(prompt=TRUE, comment=NA) #$ remove highlight and background to get standard knitr look.
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
#options(replace.assign=TRUE, width=90, prompt="R> ")
set.seed(1986)

## ----random-walk,echo=F, out.width='.6\\textwidth', message=F-------
library(randomizeR)
#R<- genSeq(crPar(10),seed=126)
R<- genSeq(crPar(10),seed=808898100)
par(cex=1.5)
plotSeq(R,plotAllSeq = T)

## ----loading--------------------------------------------------------
library("randomizeR")

## ----load-ggplot, message=FALSE, echo=FALSE-------------------------
library(ggplot2)
theme_update(text=element_text(size=24), axis.title.y=element_text(margin=margin(0,13,0,0)),  axis.title.x=element_text(margin=margin(13,0,0,0)))

## -------------------------------------------------------------------
N <- 10 
(params <- crPar(N))

## ----echo=2:3-------------------------------------------------------
set.seed(112)
params <- crPar(N)
(R <- genSeq(params))

## -------------------------------------------------------------------
getRandList(R)

## ----results=F------------------------------------------------------
saveRand(R, file="myRandList.csv")

## ----plotSeq,fig.keep='none'----------------------------------------
plotSeq(R, plotAllSeq = T)

## -------------------------------------------------------------------
(allSeqs <- getAllSeq(params))

## -------------------------------------------------------------------
N <- 50
params <- crPar(N)
(randomSeqs <- genSeq(params, r = 10000))

## ----echo=2:3-------------------------------------------------------
myPaste <- function(R) {apply(R, 1, function(x) paste(x, collapse = ""))}
p <- getProb(allSeqs)
head(data.frame(Sequences = myPaste(getRandList(allSeqs)), 
	        Probability = round(p, 6)))

## ----model, results='hide'------------------------------------------
muA <- muB <- 0
sigmaA <- sigmaB <- 1
normalEndpoint <- normEndp(mu = c(muA, muB), sigma = c(sigmaA, sigmaB))

## -------------------------------------------------------------------
(cb <- chronBias(type = "linT", theta = 1, method = "exact"))

## ----parameters-----------------------------------------------------
N <- 12
mti <- 2
bsdSeq <- getAllSeq(bsdPar(N, mti))

d <- 1.796 
sb <- selBias("CS", eta = d/4, method = "exact")
cb <- chronBias("linT", theta = 1/N, method = "exact") 
pw <- setPower(d, method = "exact") 

## ----assessment-----------------------------------------------------
(A <-assess(bsdSeq, sb, cb, pw, endp = normalEndpoint))

## ----summary-assess-------------------------------------------------
summary(A)

## ----parameters-comp------------------------------------------------
mpSeq <- getAllSeq(mpPar(N, mti))
bc <- rep(4, N/4) 
pbrSeq <- getAllSeq(pbrPar(bc))

## ----comparison-----------------------------------------------------
(C <- compare(sb, bsdSeq, mpSeq, pbrSeq, endp = normalEndpoint))

## ----plot-comparison-hide, fig.keep='none', warning=FALSE-----------
plot(C)
plot(C, y = "boxplot")

## ----plot-comparison-show-violin, echo=FALSE------------------------
plot(C)

## ----plot-comparison-show-box, echo=FALSE, warning=FALSE------------
plot(C, y='boxplot')

