## load function
source("DalbergiaWoodAnatomy_helperfunctions") # you need recode.wa() and get.diag() from here


## read example data
df <- read.delim("data/DalbergiaWoodAnatomy_exampledata.txt")


## recode from WP [ab, c, b, ...] to WP_a [0, 1] / WP_b [0, 1] / WP_c [0, 1]
# df <- dc.meta[dc.meta$Species %in% c("chlorocarpa","davidii","lemurica","urschii"),c("Species",var.wa)]
# write.table(df, file = "RawData/get.diag_exampledata.txt", row.names = T, quote = F, sep = "\t")
df.recoded <- recode.wa(df = df[,-1])


## inspect correlated variables
df.recoded$redundant # these variables are 100% correlated with other variables in the dataset

# each column and row is a variable. The number is the Pearson correlation coefficient. 1 or -1 means it's redundant.
suppressWarnings(na.omit(cor(df.recoded$df)[,df.recoded$redundant], use = "pairwise"))
# WP_b = APA_b, APP_b = APP_e, RW_a = -RW_b, RCC_a = -RCC_b, SS_a = -SS_b


## prepare dataset
View(df.recoded$df) # your recoded data
df$Species # your grouping factor. you can also use geographic regions, or anything here!
table(df$Species) # table of groups (how many of each species in this case)


## search diagnostic characters (on the basis of just the 4 example species)
res.diag <- get.diag(df = df.recoded$df, fac = df$Species) # takes a while if you have a large dataset

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
# THR threshold (1 for qualitative character states, any number for quantitative characters)
# STATE whether the GROUP has frequently a high (≥) or low (<) value for a VARIABLE
# FREQ frequency of character state VARIABLE in GROUP
# MISS percentage of missing data in VARIABLE
# TPR true positive rate (proportion of individuals of a given group showing a given character state)
# FPR false positive rate (proportion of individuals excluding a given group showing that character state)

