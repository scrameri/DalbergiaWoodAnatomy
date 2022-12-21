# DalbergiaWoodAnatomy
The R code and functions included in this repository were developed in the following article:

Ramanantsialonina, R.N., Crameri, S., Sandratriniaina, N.A., Wiemann, M.C., Hermanson, J.C., Rakouth, B., & Ravaomanalina, B.H. (2022) Comparative wood anatomy of 16 Malagasy <i>Dalbergia</i> species (Fabaceae) using multivariate techniques. IAWA J. DOI: [https://doi.org/10.1163/22941932-bja10105](https://doi.org/10.1163/22941932-bja10105)


## Folder structure
- [Ramanantisalonina_etal_2022](https://github.com/scrameri/DalbergiaWoodAnatomy/blob/main/Ramanantisalonina_etal_2022) contains the [R script](https://github.com/scrameri/DalbergiaWoodAnatomy/blob/main/Ramanantisalonina_etal_2022/Ramanantialonina_etal_2022_ANALYSES.R) for all analyses carried out in the article.

- [diagnostic_example](https://github.com/scrameri/DalbergiaWoodAnatomy/blob/main/diagnostic_example/) contains an [R script](https://github.com/scrameri/DalbergiaWoodAnatomy/blob/main/diagnostic_example/DalbergiaWoodAnatomy_example.R) demonstrating how to recode qualitative characters to dummy binary variables, combine them with quantitative characters, and search for potentially diagnostic characters to distinguish between custom groups.

## Toy dataset demonstration
The toy dataset demo is intended to show how raw data is recoded automatically to separate character state combinations, and subjected to search for characters and thresholds that are diagnostic for various groups (species or any other grouping factor), and at multiple levels (single species, group of 2 species, etc.).

### Instructions
1) Clone repository (Green button above, or click [here](https://github.com/scrameri/DalbergiaWoodAnatomy/archive/refs/heads/main.zip))
2) Make sure to install R ([Windows](https://cran.r-project.org/bin/windows/base/), [MAC](https://cran.r-project.org/bin/macosx/))
3) Optionally, also install [RStudio](https://www.rstudio.com/products/rstudio/download/)
4) Open the script [DalbergiaWoodAnatomy_example.R](https://github.com/scrameri/DalbergiaWoodAnatomy/blob/main/diagnostic_example/DalbergiaWoodAnatomy_example.R) and set the R working directory to Source File Location (RStudio: Session/Set Working Directory/To Source File Location)
5) Run it! click on first line and use cmd+Enter (Mac) or ctrl+Enter (Windows) to run it line-by-line.
