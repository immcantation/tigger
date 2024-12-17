# Download and installation

Download
-------------------------------------------------------------------------------
    
The latest stable release of TIgGER can be downloaded from 
<a href="http://cran.rstudio.com/web/packages/tigger" target="_blank">CRAN</a> or 
<a href="https://github.com/immcantation/tigger/tags" target="_blank">GitHub</a>.

Installing Released Versions
-------------------------------------------------------------------------------

The simplest way to install TIgGER is via CRAN:

```R
install.packages("tigger")
```

Downloaded source builds from GitHun may be installed in the usual way:
    
```R
install.packages("tigger_x.y.z.tar.gz", repos=NULL, type="source")
```

Building Development Versions
-------------------------------------------------------------------------------
    
To build from the [source code](https://github.com/immcantation/tigger),
first install the build dependencies:
    
```R
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown"))
```

To install the latest development code via devtools:
    
```R
library(devtools)
install_github("immcantation/tigger@master")
```

Note, using `install_github` will not build the documentation. To generate the 
documentation, clone the repository and build as normal using devtools, 
roxygen and knitr:
    
```R
library(devtools)
install_deps()
document()
build()
install()
```
