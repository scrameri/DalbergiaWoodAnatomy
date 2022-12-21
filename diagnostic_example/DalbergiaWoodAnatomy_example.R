##########################################
### Diagnostic Characters: toy example ###
##########################################

# author: sfcrameri@gmail.com, Mar 2022

## Instructions
# Set working directory to the folder where this script is located.
# You can do so by clicking on the RStudio Menu "Session/Set Working Directory/To Source File Location"

## load function
source("DalbergiaWoodAnatomy_helperfunctions.R") # you need recode.wa() and get.diag() from here

## read example data
df <- read.delim("get.diag_exampledata.txt")

## inspect
View(df)

## recode from WP [ab, c, b, ...] to WP_a [0, 1] / WP_b [0, 1] / WP_c [0, 1]
# df <- dc.meta[dc.meta$Species %in% c("chlorocarpa","davidii","lemurica","urschii"),c("Species",var.wa)]
# write.table(df, file = "RawData/get.diag_exampledata.txt", row.names = T, quote = F, sep = "\t")
df.recoded <- recode.wa(df = df[,-1]) # all except the first column (species) are wood anatomical characters

## inspect correlated variables
df.recoded$redundant # these variables are 100% correlated with other variables in the dataset

# see that WP_b = APA_b, APP_b = APP_e, RW_a = -RW_b, RCC_a = -RCC_b, SS_a = -SS_b
suppressWarnings(na.omit(cor(df.recoded$df)[,df.recoded$redundant], use = "pairwise"))
# each column and row is a variable. The number is the Pearson correlation coefficient. 1 or -1 means it's redundant.


## prepare dataset
View(df.recoded$df) # your recoded data
df$Species # your grouping factor. you can also use geographic regions, or anything here!
table(df$Species) # table of groups (how many of each species in this case)


## search diagnostic characters (on the basis of just the 4 example species)
get.diag # look at the function and Usage
res.diag <- get.diag(df = df.recoded$df, grouping = df$Species, factor.name = "Species") # takes a while if you have a large dataset

# inspect results (filter for high TPR and low FPR!)
View(res.diag) # the table contains much more than the truly diagnostic characters and thresholds, you can filter it

# filter for high true positive rate (TPR) and low false positive rate (FPR)
subset(res.diag, TPR == 1 & FPR == 0) # these are your diagnostic characters at any level. Check if it's correct!

subset(res.diag, TPR == 1 & FPR == 0 & LEVEL == 1) # these are your diagnostic characters at the level of single species. Check if it's correct!


## Explanation of output columns
# FACTOR grouping factor (e.g. species)
# NLEVELS number of grouping factor levels (e.g. number of species)
# LEVEL grouping level (equal to the length of GROUP)
# GROUP group (e.g. single species, or group of species) for which a VARIABLE is diagnostic
# VARIABLE diagnostic character (state)
# THR threshold (any number for quantitative characters, NA for qualitative character states)
# STATE whether the GROUP has frequently a high (≥) or low (<) value for a quantitative VARIABLE, or a VARIABLE is present or absent in a given GROUP
# FREQ frequency of character state VARIABLE in GROUP
# MISS percentage of missing data in VARIABLE
# TPR true positive rate (proportion of individuals of a given group showing a given character state)
# FPR false positive rate (proportion of individuals excluding a given group showing that character state)

