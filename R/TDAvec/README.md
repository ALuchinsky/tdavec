
# Overwiew

The `TDAvec` package provides implementations of several vector summaries of persistence diagrams studied in Topological Data Analysis (TDA). Each is obtained by discretizing the associated summary function computed from a persistence diagram.


# Installation

```r
# install development version from GitHub
devtools::install_github("alexey-luchinsky/TDAvec")

# install development version with vignettes/tutorials
devtools::install_github("alexey-luchinsky/TDAvec", build_vignettes = TRUE)
```

# Sample Code

Example below shows how to use the `computeVPB` function:

```r
N <- 100 
set.seed(123)
# sample N points uniformly from unit circle and add Gaussian noise
X <- TDA::circleUnif(N,r=1) + rnorm(2*N,mean = 0,sd = 0.2)

# compute a persistence diagram using the Rips filtration built on top of X
D <- TDA::ripsDiag(X,maxdimension = 1,maxscale = 2)$diagram 

# switch from the birth-death to the birth-persistence coordinates
D[,3] <- D[,3] - D[,2] 
colnames(D)[3] <- "Persistence"

# construct one-dimensional grid of scale values
ySeqH0 <- unique(quantile(D[D[,1]==0,3],probs = seq(0,1,by=0.2))) 
tau <- 0.3 # parameter in [0,1] which controls the size of blocks around each point of the diagram 
# compute VPB for homological dimension H_0
computeVPB(D,homDim = 0,xSeq=NA,ySeqH0,tau)

xSeqH1 <- unique(quantile(D[D[,1]==1,2],probs = seq(0,1,by=0.2)))
ySeqH1 <- unique(quantile(D[D[,1]==1,3],probs = seq(0,1,by=0.2)))
# compute VPB for homological dimension H_1
computeVPB(D,homDim = 1,xSeqH1,ySeqH1,tau) 
```

More information can be found in the package vignette file or help pages of the other functions.

# Citations

It you are using `TDAvec`, consider citing the article

   Chan, K. C., Islambekov, U., Luchinsky, A., & Sanders, R. (2022). A computationally efficient framework for vector representation of persistence diagrams. _Journal of Machine Learning Research_, 23, 1-33.
   
   
