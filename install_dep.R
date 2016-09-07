library(devtools)
library(versions)

tigger_version <- read.dcf("DESCRIPTION", fields = "Version")

d <- read.dcf("DESCRIPTION", fields = "Depends")
d <- sub("^\\n","",d)
depends <- strsplit(d,",\n")[[1]]
alakazam <- depends[grep("alakazam", depends)]
shazam <- depends[grep("shazam", depends)]

installDep <- function(this_pack_v, dep_pack_name, dep_pack_v) {
    required_version <- gsub(".*\\([^0-9.]*(.*)\\)$", "\\1", dep_pack_v)
    devel <- length(grep("\\.999$",required_version)) > 0
    
    cran_versions <- available.versions(dep_pack_name)
    cran_versions <- cran_versions[[dep_pack_name]]$version
    
    this_pack_devel <- length(grep("\\.999$",this_pack_v)) > 0
    
    if (!this_pack_devel & !devel & (required_version %in% cran_versions)) {
        tryCatch({install.versions(dep_pack_name,required_version)},
                 warning=function(w) { 
                     print(w)
                     message("\nInstalling from Bitbucket...\n ")
                     install_bitbucket(paste0("kleinstein/",dep_pack_name,"@default"))
                     })
    } else {
        if (!devel) { warning(paste0(required_version," not found in CRAN. Install most recent version from Bitbucket instead.")) }
        install_bitbucket(paste0("kleinstein/",dep_pack_name,"@default"))
    }
}

installDep(tigger_version, "alakazam", alakazam)
installDep(tigger_version, "shazam", shazam)

