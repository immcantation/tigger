Download
-------------------------------------------------------------------------------
    
The latest stable release of TIgGER can be downloaded from 
<a href="http://cran.rstudio.com/web/packages/tigger" target="_blank">CRAN</a> or 
<a href="https://bitbucket.org/kleinstein/tigger/downloads" target="_blank">Bitbucket</a>.

Installing Released Versions
-------------------------------------------------------------------------------

The simplest way to install TIgGER is via CRAN:

```R
install.packages("tigger")
```

Downloaded source builds from Bitbucket may be installed in the usual way:
    
```R
install.packages("tigger_x.y.z.tar.gz", repos=NULL, type="source")
```

Building Development Versions
-------------------------------------------------------------------------------
    
To build from the [source code](https://bitbucket.org/kleinstein/tigger),
first install the build dependencies:
    
```R
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown"))
```

To install the latest development code via devtools:
    
```R
library(devtools)
install_bitbucket("kleinstein/tigger@default")
```

Note, using `install_bitbucket` will not build the documentation. To generate the 
documentation, clone the repository and build as normal using devtools, 
roxygen and knitr:
    
```R
library(devtools)
install_deps()
document()
build()
install()
```
