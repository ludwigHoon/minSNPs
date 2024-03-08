# minSNPs Package
![version](https://img.shields.io/github/v/release/ludwigHoon/MinSNPs?sort=semver)
![licence](https://img.shields.io/badge/licence-MIT-blue)
![last change](https://img.shields.io/github/last-commit/ludwigHoon/MinSNPs)

***
For manual, see [here](https://github.com/ludwigHoon/minSNPs/blob/master/docs/usermanual.pdf).

This is a rebuilt version of minSNPs. Its functionality is similar to the original repository, but has cleaner code and less dependencies.

The package has a dependency that is not in CRAN, but in BioC. In order for R to install that, you need to enable BioC repository, this can be done by running `setRepositories()` and selecting both CRAN & BioC Software.

This is a simple install script for R packages:
```R
if (!require("minSNPs", quietly = TRUE)){
  if (!require("BiocParallel", quietly = TRUE)){
    if (!require("BiocManager", quietly = TRUE)){
      install.packages("BiocManager")
    }
    BiocManager::install("BiocParallel")
  }
  install.packages("minSNPs")
}
```

See [standard workflow](https://github.com/ludwigHoon/minSNPs/blob/master/minsnps_standard_workflow.R) or [standard workflow in Google colab](https://colab.research.google.com/drive/15MQFZIGzrnFy12XpUJ3345_VLqPGwgl1?usp=sharing) for quick start.