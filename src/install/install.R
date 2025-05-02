# install.packages("BiocManager")
# BiocManager::install("minfi", ask = T)
# BiocManager::install("limma", ask = T)
# BiocManager::install("IlluminaHumanMethylation450kmanifest", ask = T)
# BiocManager::install("IlluminaHumanMethylationEPICmanifest", ask = T)
# BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", ask = T)
# BiocManager::install("dplyr", ask = T)

install.packages("BiocManager", dependencies = FALSE)

packages <- c(
    "minfi",
    "limma",
    "IlluminaHumanMethylation450kmanifest",
    "IlluminaHumanMethylationEPICmanifest",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
    "dplyr"
)

BiocManager::install(packages, ask = FALSE, update = FALSE, dependencies = c("Depends", "Imports"))