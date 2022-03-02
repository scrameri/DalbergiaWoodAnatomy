############################
### ANALYZE WOOD ANATOMY ###
############################

## Libraries
need.pckg <- c("readxl","tidyverse","cluster","ape","terra","GGally","CCA","dendextend","ComplexHeatmap","plotly")
for (i in need.pckg) if (!i %in% installed.packages()) install.packages(i)

library(readxl) # read_excel
library(tidyverse) # ggplot
# library(tidymodels) # tune_grid etc.
library(cluster) # daisy
library(ape) # as.phylo
library(terra) # rast
library(GGally) # ggcorr
library(CCA) # CCA, matcor
library(dendextend) # color_branches
library(ComplexHeatmap) # BiocManager::install("ComplexHeatmap")
library(plotly) # plot_ly()
library(ggrepel)
library(ggforce)
library(rpart.plot)

## Helperfunctions
source("DalbergiaWoodAnatomy_helperfunctions.R")

## Parameters
# file paths
# input <- "RawData/20211101_Dalbergia_woodanatomy_raw.xlsx"
# input <- "RawData/20220120_Dalbergia_woodanatomy_raw_meta.xlsx"
input <- "RawData/20220223_Dalbergia_woodanatomy_raw.xlsx"
# db <- "http://raw.githubusercontent.com/scrameri/MBG/main/data/ZT_G3D_public.txt"
p.dn.meta <- "RawData/metadata_long.rds"
p.da.meta <- "RawData/metadata_atlas_long.rds"

# general
pal.species <- spectre # funky
pal.species2 <- spectre2 # funky2
scale <- TRUE
write.pdf <- TRUE

# wood anatomy variables
redundant <- c("TGR_f")
var.id <- "id"
var.excl <- c("TVD_Porous_Earlywood","TVD_Porous_Latewood","TVD_DIFFUSE")
rank <- c("GR_a"=5,"GR_b"=5,"TGR_a"=4,"TGR_b"=2,"TGR_c"=4,"TGR_d"=4,"TGR_e"=4,"TGR_f"=3,"WP_a"=1,"WP_b"=5,"WP_c"=4,"APA_a"=4,"APA_b"=4,
          "APP_a"=4,"APP_b"=4,"APP_c"=4,"APP_d"=1,"APP_e"=4,"BP_a"=4,"BP_b"=4,"BP_c"=4,"BP_d"=4,"SS_a"=1,"SS_b"=4,"SS_"=4,"RW_a"=4,"RW_"=4,
          "RW_b"=4,"VD"=1,"VL"=1,"IPD"=1,"FL"=1,"RMM"=1,"RH"=1,"NCR"=1,"NSS"=1,
          "VG_a"=1,"PP_a"=1,"VP_a"=1,"VRP_a"=1,"GTF_a"=4,"APCT_a"=1,"APCT_b"=1,"RCC_a"=4,"RCC_b"=4,"PC_a"=1,"PC_b"=1)

# ecovariables
geovars <- c("LongitudeDecimal","LatitudeDecimal")
proj.longlat <- "+init=epsg:4326"
ecovar.qual <- c("ECO_vegetationAtlas2007","ECO_surfaceLithologyAfrica","ECO_ProtectedAreas","ECO_ecoregion_Dinerstein2017")
min.frac.dummy <- 0.05 # minimum level frequency: at least this fraction of observations need to be at this ecovar.qual level for that level to be kept


################################################################################

## On the nomenclature of R objects:
# dn = New dataset (93)
# da = Atlas dataset (19)
# dd = Detienne dataset (24)
# di = InsideWood dataset (24)
# dc = New (93) + Atlas (8)
# dx = New (93) + Atlas (8) + Detienne (4) + InsideWood (4)

# dX.orig = no dummy variables
# dX = dummy variables
# dX.meta = dummy variables plus additional data (ecological, groupings, etc)


## Read New dataset
dn <- data.frame(read_excel(input, sheet = "New", na = "NA")) # rownames(dn) have ID and Taxon
rownames(dn) <- sapply(strsplit(sub("^D_|^ATLAS_D_", "", dn[,var.id]), split = "_"), "[", 2)
dn <- dn[,!names(dn) %in% c(var.id)]
dn.orig <- dn
dn.meta <- readRDS(p.dn.meta)
dn.meta[,"Species_Atlas"] <- dn.meta[,"Determination_Atlas"] <- NA # makes dn.meta conform with da.meta

## Read Atlas dataset
da <- data.frame(readxl::read_excel(input, sheet = "Atlas", na = "NA"))[,c(var.id, names(dn.orig),"Species_Atlas", "Determination_Atlas", "IAWA")]
rownames(da) <- da[,var.id]
da <- da[,!names(da) %in% c(var.id)]
da$FL <- da$NCR <- as.numeric(NA)
da.orig <- da
da.meta <- readRDS(p.da.meta)
rownames(da.meta) <- gsub("Atlas_Atlas_", "Atlas_", paste0("Atlas_", rownames(da.meta)))

# combine taxa into sensu lato
relabel <- function(x) {
  gsub("aff._hildebrandtii_N", "hildebrandtii",
       gsub("cf._purpurascens|purpurascens", "purpurascens s.l.",
            gsub("aff._normandii_Makirovana|notseen_normandii", "normandii", x)))
}
dn.meta[,"Species"] <- relabel(dn.meta[,"Species"])
da.meta[,"Species"] <- relabel(da.meta[,"Species"])

## Read Detienne / InsideWood strings
ddi <- data.frame(readxl::read_excel(input, sheet = "InsideWood", na = "NA"))

## Read Lookup table with correspondance between IAWA codes and codes in Table 1
lookup <- data.frame(readxl::read_excel(input, sheet = "IAWAConversionTable", na = c("","NA")))


message("Recode data from abc to wide format?")
scan()


## Annotate and restructure data (dummy variables)
# make dummy variables, e.g. "GR" -> c("GR_a", "GR_b")
# remove redundant variables (one out of two binary dummy variables)
var.wa <- names(dn.orig)[!names(dn.orig) %in% var.excl] 
stopifnot(length(var.wa) %in% c(25, 28)) # 25 variables in Table 1

# Ravo 93
dn.meta <- annotate.meta(dx = dn, dy = dn.meta, fac = dn.meta[,"Species"],
                         qual = ecovar.qual, min.frac.dummy = min.frac.dummy)
res.dn <- recode.wa(df = dn.orig[,var.wa])
red.dn <- c("GR_b","TGR_f","TGR_b","RW_b","SS_b") # GR_b = TGR_f = -GR_a ; TGR_b = -WP_c ; RW_b = -RW_a ; SS_b = -SS_a
# dn <- res.dn$df[,!names(res.dn$df) %in% red.dn] # 39 character states
dn <- res.dn$df[,!names(res.dn$df) %in% redundant]
sum(!names(res.dn$df) %in% red.dn) # 42 character states
dn.meta$Dataset <- "New"

# Atlas 19
# scan()
# res.da <- recode.wa(df = da[,var.wa])
# da2 <- recode.IAWA(strings = da[,"IAWA"], names = rownames(da), lookup = lookup, var.qual = var.qual, var.quant = var.quant)
# 
# names(da2$df)
# abs <- c("SS_b")
# for (i in names(da2$df)[!names(da2$df) %in% abs]) {
#   # APP_e and BP_d and RW_b and SS_b are partly assessed in da2, absent in res.da
#   if (!all(is.na(da2$df[,i])) & sum(da2$df[,i],na.rm=T) > 0) {
#     cat(i)
#     print(all.equal(as.numeric(da2$df[,i]), as.numeric(res.da$df[,i])))
#   }
# }
# 
# # SS_a for RIR2416 should have been coded as 119 absent in the Atlas
# # IPD (4-7) for RIR2418 should have been coded as 25, not 26 in the Atlas
# # RW_a for RZK7709, RBE2247, CR6512, RIR2446 and RIR2470 should have been coded as 96 in the Atlas, not only 96.1
# # RH for RIR2416 should have been coded as 102.X
# state <- "RW_b"
# dtest <- data.frame(IAWA = da2$df[,state], Ravo = res.da$df[,state],da[,c(sub("_.*$","", state),"Species_Atlas")])
# dtest$match <- ifelse(dtest[,1] == dtest[,2], 1, 0)
# dtest
# 
# write.table(da2$df, file = "RawData/20220211_Atlas_IAWA2Ravo.csv", sep = ";", dec = ".", row.names = T)

da.meta <- annotate.meta(dx = da, dy = da.meta, fac = da.meta[,"Species"],
                         qual = ecovar.qual, min.frac.dummy = min.frac.dummy)
res.da <- recode.wa(df = da[,var.wa])
names(which(apply(res.da$df, 2, function(x) var(x, na.rm = T)) == 0)) # see monomorphic
red.da <- c("TGR_f","TGR_b","RW_b") # TGR_b = -WP_c (but more missingness in TGR_b)
# da <- res.da$df[,!names(res.da$df) %in% red.da] # 36 character states
da <- res.da$df[,!names(res.da$df) %in% redundant]
sum(!names(res.da$df) %in% red.da)# 36 character states (40 with TVD)
da.meta$Dataset <- "Atlas"

# Ravo + Atlas 93+8
idx.atlas <- da.meta[,"Species"] %in% dn.meta[,"Species"]
# sort(unique(da.meta[,"Species"][idx.atlas]))
stopifnot(sum(idx.atlas) == 8) # 8 collections from 7 species
res.dc <- recode.wa(df = rbind(dn.orig, da.orig[idx.atlas,names(dn.orig)])[,var.wa])
red.dc <- c("TGR_f","TGR_b","RW_b","SS_b") # SS_a = -SS_b ; TGR_b = -WP_c (but more missingness in TGR_b)
# dc <- res.dc$df[,!names(res.dc$df) %in% red.dc] # 34 character states
dc <- res.dc$df[,!names(res.dc$df) %in% redundant]
sum(!names(res.dc$df) %in% red.dc) # 44 character states: 36 qualitative, 8 quantitative
dc.meta <- bind_rows(dn.meta[,!grepl("^TVD",names(dn.meta))], da.meta[idx.atlas,!grepl("^TVD",names(da.meta))])
names(dc.meta) <- gsub("/", "_", names(dc.meta))
dc.meta[,"Species_N"] <- factor(paste0("D. ", sub("([A-Za-z. ])_n([0-9]+)", "\\1 (n = \\2", annotate.n(dc.meta[,"Species"])), ")"))
dc.meta$Dataset <- "New"
dc.meta$Dataset[grepl("Atlas", rownames(dc.meta))] <- "Atlas"

var.qual <- names(dc)[grepl("_", names(dc)) & !grepl("^TVD", names(dc))]
var.quant <- names(dc)[!grepl("_", names(dc)) & !grepl("^TVD", names(dc))]
var.quant.tvd <- names(dc)[!grepl("_", names(dc)) | grepl("^TVD", names(dc))]

# Detienne 7
idx.det <- !grepl("\\?|_", ddi[,"Species"])
tax7 <- sort(unique(ddi[,"Species"][idx.det]))
stopifnot(sum(idx.det) == 7) # 7 species
res.det <- recode.IAWA(strings = ddi[,"WoodAnatomy_Detienne"], names = ddi[,"Species"], lookup = lookup, var.qual = var.qual, var.quant = var.quant)

rs <- array(NA, dim = c(length(res.det$list), 333), dimnames = list(names(res.det$list), 1:333))
rl <- lapply(res.det$list, function(x) {r <- as.numeric(x[rowSums(x[,2:3])==0,1])})
for (i in seq(length(rl))) {rs[i,rl[[i]]] <- 1}
names(which(colSums(rs, na.rm = T) == length(rl)))

names(which(apply(res.det$df, 2, function(x) var(x, na.rm = T)) == 0)) # see monomorphic
dd24 <- res.det$df
rownames(dd24) <- paste0("Detienne_", paste0(sapply(strsplit(sub("aff._", "aff.", rownames(dd24)), split = "_"), "[", 1)))
dd24.meta <- data.frame(Species = ifelse(sub("^Detienne_", "", rownames(dd24)) %in% tax7, sub("^Detienne_", "", rownames(dd24)), paste0(sub("^Detienne_", "", rownames(dd24)), "_MIX")), row.names = rownames(dd24))
dd24.meta$Dataset <- "Détienne"

# InsideWood 7
res.iw <- recode.IAWA(strings = ddi[,"WoodAnatomy_InsideWood"], names = ddi[,"Species"], lookup = lookup, var.qual = var.qual, var.quant = var.quant)
names(which(apply(res.iw$df, 2, function(x) var(x, na.rm = T)) == 0)) # see monomorphic
di24 <- res.iw$df
rownames(di24) <- paste0("InsideWood_", paste0(sapply(strsplit(sub("aff._", "aff.", rownames(di24)), split = "_"), "[", 1)))
di24.meta <- data.frame(Species = ifelse(sub("^InsideWood_", "", rownames(di24)) %in% tax7, sub("^InsideWood_", "", rownames(di24)), paste0(sub("^InsideWood_", "", rownames(di24)), "_MIX")), row.names = rownames(di24))
di24.meta$Dataset <- "InsideWood"

# Detienne 4
idx.iw <- ddi[,"Species"] %in% dc.meta[,"Species"]
tax4 <- sort(unique(ddi[,"Species"][idx.iw]))
stopifnot(sum(idx.iw) == 4) # 4 species
dd <- dd24[idx.iw,]
dd.meta <- dd24.meta[idx.iw,,drop=F]

# InsideWood 4
di <- di24[idx.iw,]
di.meta <- di24.meta[idx.iw,,drop=F]

# Ravo + Atlas + Detienne + InsideWood 93+8+4+4
dx <- rbind(dc[,names(dc) %in% names(dd)], dd, di)
dx.meta <- bind_rows(dc.meta, dd.meta, di.meta)
idx.x <- is.na(dx.meta[,"Genus"])
dx.meta[idx.x,"Genus"] <- "Dalbergia"
dx.meta[idx.x,"ZT_DM_SP"] <- "NA"
dx.meta[idx.x,"ZT_Dalbergia_Group"] <- rep(c("Baronii","Chlorocarpa","Monticola","Monticola"), 2)
for (i in c("Genus","ZT_DM_SP","ZT_Dalbergia_Group","ZT_Dalbergia_Supergroup","FullName")) {
  dx.meta[idx.x,i] <- dx.meta[unique(match(dx.meta[idx.x,"Species"], dx.meta[,"Species"])),i]
  dx.meta[,i] <- factor(as.character(dx.meta[,i]))
}
dx.meta[,"Species_N"] <- factor(paste0("D. ", sub("([A-Za-z. ])_n([0-9]+)", "\\1 (n = \\2", annotate.n(dx.meta[,"Species"])), ")"))
load("~/Dropbox/_ETH_Arbeit/Dalbergia/Sampling/data.rda") # load hidden coordinates
dx.meta <- cbind(dx.meta, data[(dx.meta$ID_Lab),c("LongitudeDecimal","LatitudeDecimal")])
dx.meta$Geogroup <- factor("W", levels = c("E","N","W"))
dx.meta$Geogroup[dx.meta$LatitudeDecimal > -14] <- "N"
dx.meta$Geogroup[dx.meta$ZT_Dalbergia_Group %in% c("Baronii","Maritima","Monticola")] <- "E"

# save
save(dn, dn.meta, file = "metadata_new.rda")
save(da, da.meta, file = "metadata_atlas.rda")
save(dc, dc.meta, file = "metadata_new_atlas.rda")
save(dx, dx.meta, file = "metadata_new_atlas_det_iw.rda")
var.wa47 <- unique(c(names(dn),names(da),names(dd),names(di),names(dc),names(dx)))

## Handle NA and remove monomorphic
dn.pca <- df2pca(df = dn[,!names(dn) %in% red.dn & !grepl("^TVD", names(dn))], fac = dn.meta[,"Species_N"]) # 29 / 39 kept # 20 missing values
da.pca <- df2pca(df = da[,!names(da) %in% red.da & !grepl("^TVD", names(da))], fac = da.meta[,"Species_N"]) # 24 / 37 kept # 89 missing values
dc.pca <- df2pca(df = dc[,!names(dc) %in% red.dc & !grepl("^TVD", names(dc))], fac = dc.meta[,"Species_N"]) # 31 / 41 kept # 73 missing values
dx.pca <- df2pca(df = dx, fac = dx.meta[,"Species_N"]) # 35 / 44 kept

## Subsets of variables with variation
dn.v <- dn[,apply(dn, 2, function(x) var(x, na.rm = T) > 0)]
dc.v <- dc[,apply(dc, 2, function(x) var(x, na.rm = T) > 0)]
dx.v <- dx[,apply(dx, 2, function(x) var(x, na.rm = T) > 0)]

message("Do correlation plot?")
scan()

###########################
## Correlation Structure ##
###########################

cor.thr <- 0.7

if (write.pdf) { pdf(paste0("Figures/CORRELATIONMATRIX.pdf"), width = 12, height = 12)

# New 93
print(
  ggcorr(data = dn, label = T, label_size = 3, size = 3, label_round = 1, legend.size = 9, method = c("pairwise","spearman")) +
    ggtitle(paste0("New 93 (n = ", nrow(dn), " ; p = ", ncol(dn), ")"))
  )
print(
  ggcorr(data = dn.v, label = T, label_size = 3, size = 3, label_round = 1, legend.size = 9, method = c("pairwise","spearman")) +
    ggtitle(paste0("New 93 - variant (n = ", nrow(dn.v), " ; p = ", ncol(dn.v), ")"))
)
print(
  ggcorr(data = dn.v, cor_matrix = subset.cor(dn.v, thr = cor.thr),
         label = T, label_size = 3, size = 3, label_round = 2, legend.size = 9, method = c("pairwise","spearman")) +
    ggtitle(paste0("New 93 - variant (n = ", nrow(dn.v), " ; p (>=", cor.thr, ") = ", ncol(subset.cor(dc.v, thr = cor.thr)), ")"))
)
# New 93 + Atlas 8
print(
  ggcorr(data = dc, label = T, label_size = 3, size = 3, label_round = 1, legend.size = 9, method = c("pairwise","spearman")) +
    ggtitle(paste0("New 93 + Atlas 8 (n = ", nrow(dc), " ; p = ", ncol(dc), ")"))
  )
print(
  ggcorr(data = dc.v, label = T, label_size = 3, size = 3, label_round = 1, legend.size = 9, method = c("pairwise","spearman")) +
    ggtitle(paste0("New 93 + Atlas 8 - variant (n = ", nrow(dc.v), " ; p = ", ncol(dc.v), ")"))
)
print(
  ggcorr(data = dc.v, cor_matrix = subset.cor(dc.v, thr = cor.thr),
         label = T, label_size = 3, size = 3, label_round = 2, legend.size = 9, method = c("pairwise","spearman")) +
    ggtitle(paste0("New 93 + Atlas 8 - variant (n = ", nrow(dc.v), " ; p (>=", cor.thr, ") = ", ncol(subset.cor(dc.v, thr = cor.thr)), ")"))
)
# New 93 + Atlas 8 + Detienne 4 + InsideWood 4
print(
  ggcorr(data = dx, label = T, label_size = 3, size = 3, label_round = 1, legend.size = 9, method = c("pairwise","spearman")) +
    ggtitle(paste0("New 93 + Atlas 8 + Detienne 4 + InsideWood 4 (n = ", nrow(dx), " ; p = ", ncol(dx), ")"))
)
print(
  ggcorr(data = dx.v, label = T, label_size = 3, size = 3, label_round = 1, legend.size = 9, method = c("pairwise","spearman")) +
    ggtitle(paste0("New 93 + Atlas 8 + Detienne 4 + InsideWood 4 - variant (n = ", nrow(dx.v), " ; p = ", ncol(dx.v), ")"))
)
print(
  ggcorr(data = dx.v, cor_matrix = subset.cor(dx.v, thr = cor.thr),
         label = T, label_size = 3, size = 3, label_round = 2, legend.size = 9, method = c("pairwise","spearman")) +
    ggtitle(paste0("New 93 + Atlas 8 + Detienne 4 + InsideWood 4 - variant (n = ", nrow(dx.v), " ; p (>=", cor.thr, ") = ", ncol(subset.cor(dc.v, thr = cor.thr)), ")"))
)
graphics.off() }

message("Do Ecology plot?")
scan()


###########################
## Ecology Plot ###########
###########################

### Ecology Plot
## Ecological dataset (Ravo 93 + Atlas 8)
d.eco <- cbind(dc[,!names(dc) %in% names(dc.meta)], dc.meta)
# names(d.eco) <- gsub("/", "_", names(d.eco))
var.ecoregion <- names(d.eco)[grep("^ECO_ecoregion_Dinerstein2017",names(d.eco))]
d.eco[,"ECO_ecoregion_Dinerstein2017"] <- apply(d.eco[,var.ecoregion], 1, function(x) {r <- var.ecoregion[which(x > 0)] ; if (length(r) == 0) r <- "ECO_ecoregion_Dinerstein2017_Madagascar subhumid forests" ; r})
# table(d.eco[!grepl("Atlas",rownames(d.eco)),"ECO_ecoregion_Dinerstein2017"])
d.eco$GR_a <- as.numeric(grepl("a",d.eco$GR))

# fix two false ecoregion assignment (subhumid) at ecoregion boundaries
d.eco["KAD0200","ECO_ecoregion_Dinerstein2017"] <- "ECO_ecoregion_Dinerstein2017_Madagascar humid forests"
d.eco["SH0681","ECO_ecoregion_Dinerstein2017"] <- "ECO_ecoregion_Dinerstein2017_Madagascar dry deciduous forests"
# table(d.eco[!grepl("Atlas",rownames(d.eco)),"ECO_ecoregion_Dinerstein2017"])
# table(d.eco[grepl("Atlas",rownames(d.eco)),"ECO_ecoregion_Dinerstein2017"])

# correlation
# var.eco <- names(d.eco)[grepl("^ECO_CHE_Bio01|^ECO_CHE_Bio02|^ECO_CHE_Bio03|^ECO_CHE_Bio04|^ECO_CHE_Bio12|^ECO_CHE_Bio14|^ECO_CHE_Bio15|^ECO_elevation|^ECO_aspect|^ECO_slope|^ECO_dist2", names(d.eco))]
var.eco <- names(d.eco)[grepl("^ECO_CHE_Bio01|^ECO_CHE_Bio03|^ECO_CHE_Bio04|^ECO_CHE_Bio07|^ECO_CHE_Bio12|^ECO_CHE_Bio13|^ECO_CHE_Bio14|^ECO_CHE_Bio15", names(d.eco))]

dc.eco <- subset.cor(cbind(dc[,var.wa47], d.eco[, var.eco]), thr = cor.thr)


## Plot ecology vs. growth rings
if (write.pdf) { pdf(paste0("Figures/ECOPLOT.pdf"), width = 9, height = 9)
  
  ## rasters
  # library(terra)
  # r.chelsa <- rast("~/Dropbox/_ETH_Arbeit/Dalbergia/Ecology/DataRastersVectors/Climate/CHELSA_MDG_orig.tif")
  # r.wc <- rast("~/Dropbox/_ETH_Arbeit/Dalbergia/Ecology/DataRastersVectors/Climate/WORLDCLIM_MDG_orig.tif")
  # r.alt <- rast("~/Dropbox/_ETH_Arbeit/Dalbergia/Ecology/DataRastersVectors/DIVA-GIS/MDG_alt/MDG_alt.tif")
  
  ## anova (does not fill well)
  # lm1 <- lm(asin(sqrt(dc$GR_a)) ~ (log(ECO_CHE_Bio15_prec.seasonality_kg_m2) + log(ECO_CHE_Bio04_temp.seasonality_degC) + log(ECO_CHE_Bio07_temp.range_degC))^2, data = d.eco)
  # 
  # d.eco$GR_a <- as.numeric(grepl("a", d.eco$GR))
  # lm1 <- glm(GR_a ~ (log(ECO_CHE_Bio15_prec.seasonality_kg_m2) + log(ECO_CHE_Bio04_temp.seasonality_degC) + log(ECO_CHE_Bio07_temp.range_degC))^2, family = "binomial", data = d.eco)
  # 
  # summary(lm1)
  # par(mfrow=c(2,2))
  # plot(lm1)
  # par(mfrow=c(1,1))
  
  ## growth rings and temp range
  p.eco <- cplot(d.eco,
                 x = "ECO_CHE_Bio15_prec.seasonality_kg_m2",
                 y = "ECO_CHE_Bio04_temp.seasonality_degC",
                 .hullfac = d.eco[,"ECO_ecoregion_Dinerstein2017"], palette.hull = function(n) {c("yellow","green4","orange")[1:n]},
                 hull.expand = unit(0.75, "cm"), hull.alpha = 0.2, hull.concavity = 100,
                 .fillfac = d.eco[,"Species_N"], palette.fill = pal.species,
                 .spiderfac = d.eco[,"Species_N"], palette.spider = pal.species,
                 .shapefac = factor(dc[,"GR_a"]),
                 .alphafac = interaction(factor(dc[,"GR_a"]), grepl("humid", d.eco[,"ECO_ecoregion_Dinerstein2017"])),
                 .sizefac = d.eco[,"ECO_CHE_Bio07_temp.range_degC"],
                 # .sizefac = d.eco[,"ECO_CHE_Bio03_isothermality_degC"],
                 # .sizefac = d.eco[,"ECO_elevation"],
                 points = TRUE, spider.lwd = 1.5, #alpha = 0.6, # size = 3,
                 plot = FALSE) +
    scale_shape_manual(values = c(22,21), guide = "none") + # 22 squares: no GR
    scale_alpha_manual(values = c(1,0.1,0.1,1), guide = "none")
  
  print(p.eco)
  
  ## growth rings and isothermality
  p.eco2 <- cplot(d.eco,
                  x = "ECO_CHE_Bio15_prec.seasonality_kg_m2",
                  y = "ECO_CHE_Bio03_isothermality_degC",
                  .hullfac = d.eco[,"ECO_ecoregion_Dinerstein2017"], palette.hull = function(n) {c("yellow","green4","orange")[1:n]},
                  hull.expand = unit(0.75, "cm"), hull.alpha = 0.2, hull.concavity = 100,
                  .fillfac = d.eco[,"Species_N"], palette.fill = pal.species,
                  .spiderfac = d.eco[,"Species_N"], palette.spider = pal.species,
                  .shapefac = factor(dc[,"GR_a"]),
                  .alphafac = interaction(factor(dc[,"GR_a"]), grepl("humid", d.eco[,"ECO_ecoregion_Dinerstein2017"])),
                  .sizefac = d.eco[,"ECO_CHE_Bio07_temp.range_degC"],
                  # .sizefac = d.eco[,"ECO_CHE_Bio03_isothermality_degC"],
                  # .sizefac = d.eco[,"ECO_elevation"],
                  points = TRUE, #alpha = 0.6, # size = 3,
                  plot = FALSE) +
    scale_shape_manual(values = c(22,21), guide = "none") + # 22 squares: no GR
    scale_alpha_manual(values = c(1,0.1,0.1,1), guide = "none")
  
  print(p.eco2)
  
  # Numeric
  for (v in var.quant.tvd) {
    p.ecoTS <- cplot(d.eco,
                    x = "ECO_CHE_Bio04_temp.seasonality_degC",
                    y = v,
                    .hullfac = d.eco[,"ECO_ecoregion_Dinerstein2017"], palette.hull = function(n) {c("yellow","green4","orange")[1:n]},
                    hull.expand = unit(0.25, "cm"), hull.alpha = 0.2, hull.concavity = 100,
                    .fillfac = d.eco[,"Species_N"], palette.fill = pal.species,
                    .spiderfac = d.eco[,"Species_N"], palette.spider = pal.species,
                    # .shapefac = factor(dc[,"GR_a"]),
                    # .alphafac = interaction(factor(dc[,"GR_a"]), grepl("humid", d.eco[,"ECO_ecoregion_Dinerstein2017"])),
                    # .sizefac = d.eco[,"ECO_CHE_Bio07_temp.range_degC"],
                    # .sizefac = d.eco[,"ECO_CHE_Bio03_isothermality_degC"],
                    # .sizefac = d.eco[,"ECO_elevation"],
                    points = TRUE, alpha = 0, size = 3, spider.lwd = 1.5,
                    plot = FALSE) +
      scale_shape_manual(values = c(22,21), guide = "none") + # 22 squares: no GR
      scale_alpha_manual(values = c(1,0.1,0.1,1), guide = "none") +
      ggtitle(v)
    
    p.ecoPS <- cplot(d.eco,
                    x = "ECO_CHE_Bio15_prec.seasonality_kg_m2",
                    y = v,
                    .hullfac = d.eco[,"ECO_ecoregion_Dinerstein2017"], palette.hull = function(n) {c("yellow","green4","orange")[1:n]},
                    hull.expand = unit(0.25, "cm"), hull.alpha = 0.2, hull.concavity = 100,
                    .fillfac = d.eco[,"Species_N"], palette.fill = pal.species,
                    .spiderfac = d.eco[,"Species_N"], palette.spider = pal.species,
                    # .shapefac = factor(dc[,"GR_a"]),
                    # .alphafac = interaction(factor(dc[,"GR_a"]), grepl("humid", d.eco[,"ECO_ecoregion_Dinerstein2017"])),
                    # .sizefac = d.eco[,"ECO_CHE_Bio07_temp.range_degC"],
                    # .sizefac = d.eco[,"ECO_CHE_Bio03_isothermality_degC"],
                    # .sizefac = d.eco[,"ECO_elevation"],
                    points = TRUE, alpha = 0, size = 3, spider.lwd = 1.5,
                    plot = FALSE) +
      scale_shape_manual(values = c(22,21), guide = "none") + # 22 squares: no GR
      scale_alpha_manual(values = c(1,0.1,0.1,1), guide = "none") +
      ggtitle(v)
    
    print(p.ecoTS)
    print(p.ecoPS)
    
  }

  # correlation
  '
  02 (diurnal temp range) and TGR_b and -WP_c
  03 (isothermality) and 02 and TGR_b and -WP_c and WP_b and -GR_b and GR_a and -BP_b
  14 and -04 and -01 and 
  15 and -14 and 03 and 01
  '
  print(
    ggcorr(data = cbind(dc[,var.wa47], d.eco[, var.eco]), cor_matrix = dc.eco,
           label = T, label_size = 3, size = 3, label_round = 2, legend.size = 9, method = c("pairwise","spearman")) +
      ggtitle(paste0("New 93 + Atlas 8 - variant (n = ", nrow(dc), " ; p (>=", cor.thr, ") = ", ncol(dc.eco), ")"))
  )
  graphics.off()}


message("Do CCA?")
scan()


############################
## CCA #####################
############################

# ## Load random sample of ecological variables
# load("RawData/r.sample_10000.rda") # r.sample: 10000 random samples of 19 Bioclimatic variables
# 
# for (i in ecovar.qual) {
#   if (i %in% names(r.sample)) {
#     tab.dummy <- sort(table(r.sample[,i], useNA = "always"), decreasing = TRUE)
#     
#     # keep same levels as in sample
#     # id.dummy <- names(tab.dummy[tab.dummy > sum(tab.dummy) * min.frac.dummy])
#     id.dummy <- gsub(paste0("^", i,"_"), "", names(dn.meta)[grep(i,names(dn.meta))])
#     
#     r.sample <- factor2dummy(df = r.sample, factorname = i, keep.levels = id.dummy)
#   }
# }

# var.string <- "^ECO_WOR|^ECO_elevation|^ECO_aspect|^ECO_slope|^ECO_dist|^ECO_surfaceLithology"
# var.string <- "^ECO_CHE|^ECO_surface"
# var.string <- "Bio01|Bio02|Bio03|Bio04|Bio12|Bio13|Bio15"
# var.string <- "^ECO_CHE_Bio01|^ECO_CHE_Bio02|^ECO_CHE_Bio03|^ECO_CHE_Bio04|^ECO_CHE_Bio12|^ECO_CHE_Bio13|^ECO_CHE_Bio15|^ECO_surface|^ECO_slope"
# var.string <- "^ECO_CHE_Bio01|^ECO_CHE_Bio02|^ECO_CHE_Bio03|^ECO_CHE_Bio04|^ECO_CHE_Bio12|^ECO_CHE_Bio14|^ECO_CHE_Bio15"
# var.string <- "^ECO_CHE_Bio01|^ECO_CHE_Bio02|^ECO_CHE_Bio03|^ECO_CHE_Bio04|^ECO_CHE_Bio12|^ECO_CHE_Bio14|^ECO_CHE_Bio15|^ECO_surface"
# var.y <- names(dc.meta)[grep(var.string, names(dc.meta))]

var.y <- var.eco
# var.x <- names(dc.pca)
var.x <- names(dc.pca)[grep("^GR_a|^WP_|^APA_|^RW_a|^SSa", names(dc.pca))]

## Define data
# X <- varFilter(X = dc.pca, cor.method = "spearman", cor.thr = 0.9, freqCut = 100, rank = rank) # wood anatomy Ravo 93 + Atlas 8
# X.meta <- dc.meta
# X <- varFilter(X = dn.pca[,var.x], cor.method = "spearman", cor.thr = 0.9, freqCut = 100, rank = rank) # wood anatomy Ravo 93
# X <- varFilter(X = dc.pca[,var.x], cor.method = "spearman", cor.thr = 0.7, freqCut = 100, rank = rank) # wood anatomy Ravo 93
X <- dc.pca[,var.x]
X.meta <- dc.meta


Y <- dc.meta[,var.y] # environment Ravo 93 + Atlas 8
# Y <- dn.meta[,var.y] # environment Ravo 93 + Atlas 8

# ecorank <- data.frame(ecovar = names(Yvar), rank = c(1,2,3,1,4,4,4,4,4,4,4,1,4,4,1,4,4,4,4,1,1,1,1,1,1,1,1,1,1,1,1,1))
# Y <- varFilter(X = Y, cor.method = "spearman", cor.thr = 0.7, freqCut = 100, rank = NULL) # Bio09, 10, 11, 16, 15, 17, 18, 19 correlated

# Yvar <- r.sample[,!names(r.sample) %in% geovars]
# Yvar <- r.sample[,var.y]
# ecorank <- data.frame(ecovar = names(Yvar), rank = c(1,2,3,1,4,4,4,4,4,4,4,1,4,4,1,4,4,4,4,1,1,1,1,1,1,1,1,1,1,1,1,1))
# ranky <- ecorank$rank ; names(ranky) <- ecorank$ecovar
# Yvar <- varFilter(X = Yvar, cor.method = "spearman", cor.thr = 0.8, freqCut = 1E06, rank = ranky)
# range(cor(Yvar, method = "spearman")[lower.tri(cor(Yvar))])
# Y <- Y[,names(Yvar)[names(Yvar) %in% names(Y)]]
# Y <- Y[,names(Yvar)[grep("^Bio", names(Yvar))]]

## CCA::cc function
# https://stats.idre.ucla.edu/r/dae/canonical-correlation-analysis/
# pal.eco <- function(n) c("yellow","green4","darkgreen","orange")[1:n]

## find variables leading to non-positive definite error
# for (probs in names(X)) {
#   cc1 <- try(cc(scale(X), scale(Y)), silent = TRUE)
#   if (inherits(cc1, "try-error")) {
#     warning("removing variable(s) ", paste(probs, collapse = ","), " due to non positive definite error")
#     cc1 <- try(cc(scale(X[,!colnames(X) %in% probs]), scale(Y)[,!colnames(Y) %in% probs]), silent = TRUE)
#   }
#   if (inherits(cc1, "try-error")) {print("This did not help")} else {print(paste("exclude", probs))}
# }
# probs <- "TGR_c"
cc1 <- try(cc(scale(X), scale(Y)), silent = TRUE)
if (inherits(cc1, "try-error")) {
      warning("removing variable(s) ", paste(probs, collapse = ","), " due to non positive definite error")
      cc1 <- try(cc(scale(X[,!colnames(X) %in% probs]), scale(Y)[,!colnames(Y) %in% probs]), silent = TRUE)
}

cc1$cor
cc1$xcoef[,c(1,2)]
cc1$xcoef[which.max(abs(cc1$xcoef[,1])),1]
cc1$xcoef[which.max(abs(cc1$xcoef[,2])),2]

cc1$ycoef[,c(1,2)]

# plot(dn$SS_a, dn.meta$ECO_CHE_Bio04_temp.seasonality_degC)

# plt.cc(cc1)
# cc1$scores$corr.Y.xscores

if (write.pdf) { pdf(paste0("Figures/CCA_Dalbergia_", paste(nrow(X), ncol(X), ncol(Y), sep = "_"), ".pdf"),
    width = 9, height = 9)

flipx = -1
flipy = -1

p.cca1 <- cplot(cc1$scores$xscores,
      flipx = flipx, flipy = flipy,
      .fillfac = X.meta[,"Species_N"], palette.fill = pal.species,
      .hullfac = X.meta[,"Species_N"], palette.hull = pal.species, hull.labels = T,
      hull.labels.cex = 8,
      plot = F, size = 4, alpha = 0.6,
) +
  labs(x = paste0("CCA Dimension 1 (", round(100*cc1$cor[1]), "%)"),
       y = paste0("CCA Dimension 2 (", round(100*cc1$cor[2]), "%)"), fill = "")

p.cca2 <- cplot(cc1$scores$xscores, x = 3, y = 4,
      flipx = 1, flipy = 1,
      .fillfac = X.meta[,"Species_N"], palette.fill = pal.species,
      .hullfac = X.meta[,"Species_N"], palette.hull = pal.species, hull.labels = T,
      hull.labels.cex = 8,
      plot = F, size = 4, alpha = 0.6,
) +
  labs(x = paste0("CCA Dimension 3 (", round(100*cc1$cor[3]), "%)"),
       y = paste0("CCA Dimension 4 (", round(100*cc1$cor[4]), "%)"), fill = "")

print(p.cca1)
print(p.cca2)

cory <- data.frame(cc1$scores$corr.Y.yscores)
corx <- data.frame(cc1$scores$corr.X.yscores)

p.cca3 <- ggplot(cory) +
  geom_segment(aes(x = 0, xend = flipx*X1, y = 0, yend = flipy*X2), arrow = arrow()) +
  geom_label(aes(x = flipx*X1, y = flipy*X2, label = rownames(cory))) +
  lims(x = c(-1,1), y = c(-1,1)) +
  labs(x = paste0("CCA Dimension 1 (", round(100*cc1$cor[1]), "%)"),
       y = paste0("CCA Dimension 2 (", round(100*cc1$cor[2]), "%)")) +
  theme_bw() +
  ggtitle("Correlation between Canonical Variables and Bioclimate")

# ggplot(data.frame(cc1$xcoef)) +
#   geom_segment(aes(x = 0, xend = flipx*X1, y = 0, yend = flipy*X2), arrow = arrow()) +
#   geom_label(aes(x = flipx*X1, y = flipy*X2, label = rownames(cc1$xcoef))) +
#   #lims(x = c(-1,1), y = c(-1,1)) +
#   labs(x = paste0("CCA Dimension 1 (", round(100*cc1$cor[1]), "%)"),
#        y = paste0("CCA Dimension 2 (", round(100*cc1$cor[2]), "%)")) +
#   theme_bw() +
#   ggtitle("Correlation between Canonical Variables and Bioclimate")

p.cca4 <- ggplot(corx) +
  geom_segment(aes(x = 0, xend = flipx*X1, y = 0, yend = flipy*X2), arrow = arrow()) +
  geom_label(aes(x = flipx*X1, y = flipy*X2, label = rownames(corx))) +
  lims(x = c(-1,1), y = c(-1,1)) +
  labs(x = paste0("CCA Dimension 1 (", round(100*cc1$cor[1]), "%)"),
       y = paste0("CCA Dimension 2 (", round(100*cc1$cor[2]), "%)")) +
  theme_bw() +
  ggtitle("Correlation between Canonical Variables and Wood Anatomy")

# ggplot(data.frame(cc1$ycoef)) +
#   geom_segment(aes(x = 0, xend = flipx*X1, y = 0, yend = flipy*X2), arrow = arrow()) +
#   geom_label(aes(x = flipx*X1, y = flipy*X2, label = rownames(cc1$ycoef))) +
#   #lims(x = c(-1,1), y = c(-1,1)) +
#   labs(x = paste0("CCA Dimension 1 (", round(100*cc1$cor[1]), "%)"),
#        y = paste0("CCA Dimension 2 (", round(100*cc1$cor[2]), "%)")) +
#   theme_bw() +
#   ggtitle("Correlation between Canonical Variables and Bioclimate")

print(p.cca3)
print(p.cca4)

graphics.off() }


### 3D CCA plot
# d <- data.frame(cc1$scores$xscores)
# names(d) <- paste0("CCA_Dimension_", seq(ncol(d)))
# d[,"Species"] <- X.meta[,"Species"]
# 
# fig <- plot_ly(d, x = ~CCA_Dimension_1, y = ~CCA_Dimension_2, z = ~CCA_Dimension_3,
#                mode = "markers",
#                color = d[,"Species"], colors = pal.species(nlevels(d[,"Species"])),
#                type="scatter3d")
# 
# fig


message("Do Parallel Coordinates Plot?")
scan()

################################
### Parallel Coordinate Plot ###
################################

tvd <- TRUE

# reshape to long format and group by Collection
if (tvd) {
  quant <- names(dc)[!grepl("_", names(dc)) | grepl("^TVD_", names(dc))] # w/ TVD
} else {
  quant <- names(dc)[!grepl("_", names(dc))] # w/o TVD
}
stopifnot(length(quant) %in% c(8,11)) # 8 quantitative variables

# Select dataset
# PC <- cbind(dc.meta[,c("ID_Lab","Species_N","Dataset")],dc.pca[,quant]) # New (93) + Atlas (8) # missing values imputed
PC <- cbind(dc.meta[,c("ID_Lab","Species_N","Dataset")],dc[,quant]) # New (93) + Atlas (8) with missing values

# scale to between 0 and 100
to <- c(0,100)
for (i in quant) PC[,paste0(i, "_orig")] <- PC[,i]
PC[,quant] <- apply(PC[,quant], 2, scales::rescale, to = to)

# reshape
value_orig <- (PC[,c("ID_Lab","Species_N","Dataset",paste0(quant, "_orig"))] %>% 
    reshape2::melt(id = c("ID_Lab", "Species_N","Dataset")) %>% 
    dplyr::select(value))[,1]

PC <- PC[,c("ID_Lab","Species_N","Dataset",quant)] %>% 
  reshape2::melt(id = c("ID_Lab", "Species_N","Dataset")) # %>% 
  # group_by(ID_Lab)

# # group species by region
# PC$group <- factor("E", levels = c("E","N","W"))
# PC$group[grep("abrahamii|urschii|pseudobaronii|hirtipetala|obtusa", PC[,"Species_N"])] <- "N"
# PC$group[grep("greveana|purpurascens|chlorocarpa|davidii|lemurica", PC[,"Species_N"])] <- "W"

# order variables (x axis)
if (tvd) {
  PC$variable <- factor(PC$variable, levels = c("IPD","VD","VL","RH","RMM","NSS","NCR","FL","TVD_RingPorousEW","TVD_RingPorousLW","TVD_DiffusePorous"))
} else {
  PC$variable <- factor(PC$variable, levels = c("IPD","VD","VL","RH","RMM","NSS","NCR","FL"))
}

# aggregate
PCn <- subset(PC, Dataset == "New")
PCa <- subset(PC, Dataset == "Atlas")
# PCn.mean <- aggregate(PCn$value, by = list(PCn$Dataset, PCn[,"Species_N"], PCn$group, PCn$variable), FUN = mean, na.rm = TRUE)
# PCa.mean <- aggregate(PCa$value, by = list(PCa$Dataset, PCa[,"Species_N"], PCa$group, PCa$variable), FUN = mean, na.rm = TRUE)
# names(PCn.mean) <- names(PCa.mean) <- c("Dataset","Species_N","group","variable","value")

PCn.mean <- aggregate(PCn$value, by = list(PCn$Dataset, PCn[,"Species_N"], PCn$variable), FUN = mean, na.rm = TRUE)
PCa.mean <- aggregate(PCa$value, by = list(PCa$Dataset, PCa[,"Species_N"], PCa$variable), FUN = mean, na.rm = TRUE)
names(PCn.mean) <- names(PCa.mean) <- c("Dataset","Species_N","variable","value")

# ranges
quant.ranges <- sapply(levels(PCn$Species_N), FUN = function(x) {
  y <- subset(PCn, Species_N %in% x)
  tapply(y$value, y$variable, FUN = function(z) max(z,na.rm=T)-min(z,na.rm=T))
})
# apply(quant.ranges, MARGIN = 1, FUN = function(x) {ifelse(x >= 80, 999, 0)})
# apply(quant.ranges, MARGIN = 1, FUN = function(x) {ifelse(x >= 90, 999, 0)})

if (write.pdf) { pdf(paste0("Figures/PCP_Dalbergia_", paste(nrow(PC), nlevels(PC$variable), sep = "_"), ".pdf"),
    width = 9, height = 9)

p.pcp <- ggplot(PCn, aes(x = variable, y = value)) +
  
  # New
  geom_line(aes(group = ID_Lab, color = Species_N), size = 0.2) +
  geom_line(data = PCn.mean, aes(group = Species_N, color = Species_N), size = 1.5) +
  geom_point(aes(group = ID_Lab, color = Species_N), size = 0.5) +
  
  # Atlas
  geom_line(data = PCa, aes(group = ID_Lab), color = "gray30", size = 0.2) +
  geom_line(data = PCa.mean, aes(group = Species_N), color = "gray30", size = 0.75) +
  geom_point(data = PCa, aes(group = ID_Lab), color = "gray30", size = 0.2) +
  
  # Legend
  # geom_hline(aes(yintercept = 20), size = 0.5) +
  # geom_hline(aes(yintercept = 80), size = 0.5) +
  # geom_hline(aes(yintercept = 5), size = 0.5) +
  # geom_hline(aes(yintercept = 95), size = 0.5) +
  scale_color_manual(values = pal.species(nlevels(PCn[,"Species_N"])), guide = "none") +
  # facet_wrap(~group) +
  facet_wrap(~Species_N) +
  labs(x = "Quantitative character", y = "Relative value") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

print(p.pcp)

graphics.off() }


message("Do HEATMAP?")
scan()


###############
### HEATMAP ###
###############

# Select dataset
# HM <- dx
HM <- dx.v[,var.qual[var.qual %in% names(dx.v)]] # only qualitative variables with variation
HM.meta <- dx.meta

# # relevel quantitative
# quant <- names(HM)[!grepl("_", names(HM))]
# HM[,quant] <- apply(HM[,quant], MARGIN = 2, FUN = scales::rescale, to = c(-2,2))

cluster.rows <- F
cluster.cols <- F

# # relevel according to 3 groups (E, N, W)
# lev <- HM.meta[,"Species_N"]
# level.order <- c(2,3,9,10,12,15, 1,7,11,13,16, 4,5,6,8,14)
# lev.col <- pal.species(nlevels(lev))[level.order]

# # relevel for dc (E, N, W)
# lev <- HM.meta[,"Species_N"]
# level.order <- c(4,5,17,18,20,23, 1,14,19,21,26, 10,11,12,15,22, 2,3,6,7,8,9,13,16,24,25)
# lev.col <- c(pal.species(16)[c(2,3,9,10,12,15, 1,7,11,13,16, 4,5,6,8,14)], rep("gray45",10))
# lev <- factor(lev, levels = levels(lev)[level.order])

# # relevel for all (E, N, W)
# lev.iw <- HM.meta[,"Species_N"]
# lev <- factor(c(as.character(dc.meta[,"Species"]), lev.iw, lev.iw), levels = unique(c(unique(sort(dc.meta[,"Species"]))[c(4,5,17,18,20,23, 1,14,19,21,26, 10,11,12,15,22, 2,3,6,7,8,9,13,16,24,25)], lev.iw)))
# level.order <- c(4,5,17,18,20,23, 1,14,19,21,26, 10,11,12,15,22, 2,3,6,7,8,9,13,16,24,25, 27:46)
# lev.col <- c(c(pal.species(16)[c(2,3,9,10,12,15, 1,7,11,13,16, 4,5,6,8,14)], rep("gray45",10)), rep("gray45", 30))

lev <- factor(HM.meta[,"Species"])
nlev <- nlevels(lev)
lev.col <- pal.species(nlev)

la <- rowAnnotation(foo = anno_block(gp = gpar(fill = lev.col),
                                     labels = levels(lev), 
                                     labels_gp = gpar(col = "white", fontsize = 0)))
clev <- factor(sapply(strsplit(names(HM), split = "_"), "[", 1), levels = unique(sapply(strsplit(names(HM), split = "_"), "[", 1)))
ta <- columnAnnotation(foo = anno_block(gp = gpar(),
                                        labels = levels(clev), 
                                        labels_gp = gpar(col = "black", fontsize = 6)))

# cluster columns and rows
if (cluster.rows) {
  row_dend <- hclust(dist(HM, method = "euclidean"), method = "ward.D2")
  k.row <- 5
}
if (cluster.cols) {
  col_dend <- hclust(dist(t(HM), method = "euclidean"), method = "ward.D2")
  k.col <- 2
}

# col.order <- c(col.qual, col.quant)
col.order <- names(HM)

if (write.pdf) { pdf(paste0("Figures/HEATMAP_Dalbergia_", paste(nrow(HM), ncol(HM), sep = "_"), ".pdf"),
    width = 9, height = 9)

p.hm <- Heatmap(
  matrix = as.matrix(HM[,col.order]),
  # heatmap color
  col = circlize::colorRamp2(c(-2, -1, 0, 1, 2), c(scales::muted("darkred"), scales::muted("red", l = 50), "white", scales::muted("blue", l = 50), scales::muted("darkblue"))),
  
  # row splitting / clustering
  row_split = if (cluster.rows) NULL else lev,
  left_annotation = if (cluster.rows) NULL else la,
  cluster_rows = if (cluster.rows) dendextend::color_branches(row_dend, k = k.row) else FALSE,
  # cluster_rows = FALSE, row_order = sapply(dd.heatmap$ID_Lab, FUN = function(x) which(dd.hclust$ID_Lab == x))
  # row_order = if (cluster.rows) seq(nrow(d.heatmap)) else sapply(dd.hclust[,var.label], FUN = function(x) which(dd.heatmap[,var.label] == x)),
  row_order = seq(nrow(HM)),
  
  # column splitting / clustering
  column_split = if (cluster.cols) NULL else clev,
  top_annotation = if (cluster.cols) NULL else ta,
  cluster_columns = if (cluster.cols) dendextend::color_branches(col_dend, k = k.col) else FALSE,
  
  # columns
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  show_column_dend = FALSE,
  
  # rows
  row_labels = rownames(HM),
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 5, col = lev.col[lev]), #pal.species(nlev)[level.order][lev]),
  row_names_rot = 0,
  show_row_dend = FALSE,
  
  # row and column titles
  row_title_side = "left",
  row_title_gp = gpar(fontsize = 10),
  row_title_rot = 0, # 45 is not supported in 2.6.2
  
  column_title_side = "top",
  column_title_gp = gpar(fontsize = 0),
  column_title_rot = 90,
  
  # general
  rect_gp = gpar(col = "gray90", lwd = 0.1),
  border = "gray45",
  heatmap_legend_param = list(title = "Scaled\nValue"),
  show_heatmap_legend = TRUE,
  name = "Score"
)

print(p.hm)

graphics.off() }


message("Do PCA?")
scan()


################################
## PCA #########################
################################

# # ecological
# PCA <- dc.meta[,grepl("^ECO_",names(dc.meta)) & !grepl("^ECO_WOR_|^ECO_imputed$|^ECO_ProtectedAreas_|^ECO_ecoregion_Dinerstein2017_|^ECO_vegetationAtlas", names(dc.meta))]
# PCA$ECO_forestCover2017 <- as.numeric(as.character(PCA$ECO_forestCover2017))
# PCA.meta <- dc.meta
# 
# names(PCA.meta) <- gsub("/", "_", names(PCA.meta))
# for (ecovar in c("ECO_CHE_Bio01_temp.annual_degC","ECO_WOR_Bio01_temp.annual_degC",
#                  "ECO_CHE_Bio12_prec.annual_kg_m2","ECO_WOR_Bio12_prec.annual_kg_m2",
#                  "ECO_CHE_Bio15_prec.seasonality_kg_m2","ECO_WOR_Bio15_prec.seasonality_kg_m2")) {
#   p <- ggplot(PCA.meta, aes_string(x = ecovar, y = "Species_N")) +
#     geom_boxplot() +
#     theme_bw()
#   
#   print(p)
#   scan()
#   
# }

# wood anatomy
# PCA <- dc.pca # NCR missing in RIR3571, 67 missing in Atlas specimens (TGR, APA, APP, BP, RW, FL, NCR)
# PCA.meta <- dc.meta
PCA <- dn.pca # NCR missing in RIR3571
PCA.meta <- dn.meta

# annotate data
PCA.meta[,"Species_COL"] <- PCA.meta[,"Species"]
PCA.meta[PCA.meta[,"Dataset"] %in% "Atlas","Species_COL"] <- "Atlas"
PCA.meta[,"Species_COL"] <- factor(PCA.meta[,"Species_COL"], levels = c(unique(PCA.meta[,"Species"]), "Atlas"))

# PCA
d.pca <- prcomp(PCA, center = TRUE, scale. = scale)

# primarily separated
apply(d.pca$rotation[,1:4], 2, function(x) {id <- which.max(abs(x)) ; data.frame(var = names(x)[id], score = x[id])})

# ANOVA of PC1 and 2
PCA.meta <- factor2dummy(PCA.meta, factorname = "Species")
PCA.meta$GRdistinct_narrowBP <- interaction(PCA$GR_a, PCA$BP_a)
levels(PCA.meta$GRdistinct_narrowBP) <- c("none","GR_a","BP_a","GR_a+BP_a")

d.aov <- data.frame(d.pca$x, PCA.meta)
aov1 <- aov(PC1 ~ Species_N, data = d.aov)
summary(aov1)

# library(multcomp)
# summary(glht(aov1, linfct = mcp(Species_N = "Tukey")))

d.tukey <- data.frame(TukeyHSD(aov1)[[1]])
d.tukey$comparison <- factor(rownames(d.tukey), levels = rev(rownames(d.tukey)))

pdf("Figures/TukeyPairwiseComparison.pdf", 12, 12)

par(mfrow=c(2,2))
plot(aov1)
par(mfrow=c(1,1))

ggplot(d.tukey) +
  geom_vline(aes(xintercept = 0), color = "gray25", size = 1) +
  geom_point(aes(x = diff, y = comparison)) +
  geom_segment(aes(x = lwr, xend = upr, y = comparison, yend = comparison)) +
  labs(x = "Difference in the observed means of first principal wood anatomical component", y = "Pairwise Comparison") +
  theme_bw()

graphics.off()



if (write.pdf) { pdf(paste0("Figures/PCA_Dalbergia_", paste(nrow(PCA), ncol(PCA), sep = "_"), ".pdf"),
    width = 9, height = 9)

flipx <- -1
flipy <- -1

# Figure 6
d.pc <- reshape2::melt(d.aov[,c("Species_N","PC1","PC2")], id = "Species_N")
ggplot(d.pc, aes(x = Species_N, y = value, alpha = variable, fill = Species_N)) +
  geom_boxplot() +
  coord_flip() +
  scale_fill_manual(values = pal.species(nlevels(PCA.meta$Species_N)), guide = "none") +
  scale_alpha_manual(values = seq(1,0.98, length.out = nlevels(d.pc$variable)), guide = "none") +
  labs(x = "", y = paste0("PC 1 (", round(100*(d.pca$sdev^2/sum(d.pca$sdev^2))[1]), "%)")) +
  theme_bw()


p.pca1 <- cplot(d.pca, x = 1, y = 2,
      flipx = flipx, flipy = flipy,
      .shapefac = factor(PCA[,"GR_a"]),
      .sizefac = factor(PCA[,"BP_a"]),
      .fillfac = PCA.meta[,"Species_COL"], palette.fill = pal.species2,
      .hullfac = PCA.meta[,"Species_N"], palette.hull = pal.species, hull.labels = T, hull.labels.cex = 12,
      alpha = 0.6, plot = FALSE) + # size = 4, 
  scale_shape_manual(values = c(22,21)) +
  scale_size_manual(values = c(2,4))

p.pca2 <- cplot(d.pca, x = 3, y = 4,
      flipx = flipx, flipy = flipy,
      .shapefac = factor(PCA[,"APA_b"]),
      .sizefac = factor(PCA[,"APP_e"]),
      .fillfac = PCA.meta[,"Species_COL"], palette.fill = pal.species2,
      .hullfac = PCA.meta[,"Species_N"], palette.hull = pal.species, hull.labels = T, hull.labels.cex = 12,
      alpha = 0.6, plot = FALSE) + # size = 4, 
  scale_shape_manual(values = c(22,21)) +
  scale_size_manual(values = c(2,4))

print(p.pca1)
print(p.pca2)

d.pca$var <- round(100*(d.pca$sdev^2 / sum(d.pca$sdev^2)), 2)

barplot(d.pca$var[1:10], main = "", col = c("gray45","gray45",rep(NA,8)))
barplot(d.pca$var[1:10], main = "", col = c(rep(NA,2),"gray45","gray45",rep(NA,6)))

p.pca3 <- ggplot(data.frame(d.pca$rotation)) +
  geom_segment(aes(x = 0, xend = flipx*PC1, y = 0, yend = flipy*PC2), arrow = arrow()) +
  geom_label(aes(x = flipx*PC1, y = flipy*PC2, label = rownames(d.pca$rotation))) +
  lims(x = c(-.5,.5), y = c(-.5,.5)) +
  labs(x = paste0("PC 1 (", d.pca$var[1], "%)"), y = paste0("PC 2 (", d.pca$var[2], "%)")) +
  theme_bw() +
  ggtitle("Loadings of Wood Anatomy to PCA Axes")

p.pca4 <- ggplot(data.frame(d.pca$rotation)) +
  geom_segment(aes(x = 0, xend = 1*PC3, y = 0, yend = 1*PC4), arrow = arrow()) +
  geom_label(aes(x = 1*PC3, y = 1*PC4, label = rownames(d.pca$rotation))) +
  lims(x = c(-.5,.5), y = c(-.5,.5)) +
  labs(x = paste0("PC 3 (", d.pca$var[3], "%)"), y = paste0("PC 4 (", d.pca$var[4], "%)")) +
  theme_bw() +
  ggtitle("Loadings of Wood Anatomy to PCA Axes")

print(p.pca3)
print(p.pca4)

graphics.off() }


message("Do HCLUST?")
scan()

####################
## HCLUST ##########
####################

## Distance (Gower)
# HC <- dc[,names(dc.pca)] # New (93) + Atlas (8) with missing data
# HC.meta <- dc.meta
HC <- dn.pca # same as PCA
HC.meta <- dn.meta
rownames(HC) <- rownames(HC.meta) <- paste0("D. ", HC.meta[,"Species"], " ", rownames(HC.meta))

HC <- HC[,!apply(HC, 2, function(x) {anyNA(x)})] # remove variables with missingness

binary <- names(HC)[grep("_", names(HC))]
nominal <- names(HC)[grep("_", names(HC), invert = TRUE)]
d.daisy <- daisy(HC, metric = "gower", stand = TRUE,
                 type = list(asymm = binary, logratio = nominal))

## Hclust
hc <- hclust(d.daisy, method = "ward.D2")
# hc <- hclust(d.daisy, method = "complete")

if (write.pdf) { pdf(paste0("Figures/HCLUST_NJ_", nrow(HC), "_", ncol(HC), ".pdf"),
    width = 9, height = 9)

# http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning
cplot.phylo(as.phylo(hc), df = HC.meta, type = "phylogram", palette.edge = pal.species, palette.tip = pal.species,
            label.offset = 0.01, tip.color = "Species_N", tip.cex = 0.4,
            legend.tip = F)

## NJ
# d.nj <- nj(X = d.daisy)
d.nj <- nj(X = dist(scale(HC), method = "euclidean"))

cplot.phylo(d.nj, df = HC.meta,
            edge.color = "Species_N", edge.lwd = 2,
            tip.color = "Species_N",
            tip.label = T, label.offset = 0.2, lab4ut = "axial",
            legend.edge.cex = 1)

graphics.off() }



message("do rpart?")
scan()

###########################
## RPART ##################
###########################

# generate edge cases for Detienne datset

# rpart
RP <- dc # New (93) + Atlas (8), all 44 variables, missing values included
RP.meta <- dc.meta

# all
idx.rp <- rownames(RP)

var.rp <- names(RP)

# # rpart subset
# sort(rank)
# var.rp <- c("WP_a","APP_d","SS_a","VD","VL","IPD","FL","RMM","RH","NCR","NSS","VG_a",
#             "PP_a","VP_a","VRP_a","APCT_a","APCT_b","PC_a","PC_b","TGR_b") # rank easy
# idx.rp <- !RP.meta$Species %in% c("chlorocarpa","davidii","lemurica","urschii")
# idx.rp <- RP.meta$Species %in% c("abrahamii","pseudobaronii","monticola","greveana","purpurascens s.l.","orientalis")
# idx.rp <- RP.meta$Species %in% c("abrahamii","monticola","greveana","purpurascens s.l.","orientalis")
# idx.rp <- RP.meta$Species %in% c("orientalis","razakamalalae","baronii","bathiei","hirtipetala","obtusa")
# idx.rp <- RP.meta$Species %in% c("razakamalalae","baronii","bathiei","hirtipetala","obtusa")
# idx.rp <- RP.meta$Species %in% c("orientalis","baronii","bathiei","hirtipetala")
# idx.rp <- RP.meta$Species %in% c("abrahamii","greveana","purpurascens s.l.")
# idx.rp <- RP.meta$Species %in% c("greveana","purpurascens s.l.")
# idx.rp <- !RP.meta$Species %in% c("chlorocarpa","davidii","lemurica","urschii","normandii")
# idx.rp <- RP.meta$Species %in% c("monticola","orientalis")
# idx.rp <- RP.meta$Species %in% c("razakamalalae","normandii")
# idx.rp <- RP.meta$Species %in% c("razakamalalae","obtusa")
# idx.rp <- RP.meta$Species %in% c("baronii","hirtipetala")
# idx.rp <- RP.meta$Species %in% c("orientalis","bathiei")

# select dataset
idx.rp <- rownames(RP.meta)
var.rp <- colnames(RP)
d.rpart <- data.frame(fac = factor(RP.meta[idx.rp,"Species"]), RP[idx.rp,var.rp]) # all

rp <- do.rpart(df = d.rpart, 
               fac = "fac", nmin = min(table(d.rpart$fac)),  # 5
               control = list(minsplit = 8, # min. nb. obs. for split
                              minbucket = 1, # min. nb. obs. in leaf
                              cp = 0.01,
                              maxcompete = 4,
                              maxsurrogate = 5,
                              usesurrogate = 2,
                              surrogatestyle = 0,
                              maxdepth = 30,
                              xval = 10))

rp$accuracy # almost perfect classification
rp$pred <- predict(rp$rp.final, RP[idx.rp,var.rp], type = "class")
rp$nwrong <- nrow(d.rpart) - sum(diag(table(rp$pred, d.rpart$fac)))

if (write.pdf) { pdf(paste0("Figures/RPART_Dalbergia_", nrow(RP), "_", ncol(RP), ".pdf"),
    width = 7, height = 7)

prp(rp$rp.final, type = 3, extra = 0, # or use extra = 1
    roundint = FALSE,
    box.col = transp(pal.species(nlevels(d.rpart$fac))[rp$rp.final$frame$yval], .6),  
    fallen.leaves = F, clip.right.labs = T, varlen = 0, space = 0, cex = 0.5, 
    gap = 8, ygap = 3, xflip = F, boxes.include.gap = T,
    under = T, branch.lty = 2, split.fun = split.fun)

title(sub = paste0("Accuracy of this tree on all samples used for training: ", round(100*rp$accuracy, 2), "% (wrongly classified: ", rp$nwrong, ")"))

graphics.off() }


message("Do DIAGNOSTIC analysis?")
scan()

##################
### DIAGNOSTIC ###
##################

DG <- dx
DG.meta <- dx.meta

# select some ecological classes
DG.meta[,"ECO_CHE_Bio01_temp.annual_>23_degC"] <- as.numeric(DG.meta$ECO_CHE_Bio01_temp.annual_degC > 23)
DG.meta[,"ECO_CHE_Bio15_prec.seasonality>100_kg_m2"] <- as.numeric(DG.meta$ECO_CHE_Bio15_prec.seasonality_kg_m2 > 100)
DG.meta[,"ECO_CHE_Bio12_prec.annual_kg>2500_kg_m2"] <- as.numeric(DG.meta$ECO_CHE_Bio12_prec.annual_kg_m2 > 2500)

table(DG.meta[,"Species"], DG.meta[,"Geogroup"])
# sp <- SpatialPointsDataFrame(coords = DG.meta[!is.na(DG.meta$LongitudeDecimal),c("LongitudeDecimal","LatitudeDecimal")],
#                              data = DG.meta[!is.na(DG.meta$LongitudeDecimal),], proj4string = CRS(proj.longlat))
# get.leaflet(sp, classfac = "Species", layerfac = "Geogroup")
# get.leaflet(sp, classfac = "Species", layerfac = "Species")

diag.spec <- get.diag(df = DG, fac = factor(DG.meta$Species)) # 421 (level 8) or 736 (complete)
diag.geo <- get.diag(df = DG, fac = factor(DG.meta$Geogroup)) # 110 (complete)
diag.sg <- get.diag(df = DG, fac = factor(DG.meta$ZT_Dalbergia_Supergroup)) # 63 (complete)
diag.g <- get.diag(df = DG, fac = factor(DG.meta$ZT_Dalbergia_Group)) # 320 (complete)

diag.temp <- get.diag(df = DG, fac = factor(DG.meta[,"ECO_CHE_Bio15_prec.seasonality>100_kg_m2"]))
diag.prec.seasonality <- get.diag(df = DG, fac = factor(DG.meta[,"ECO_CHE_Bio15_prec.seasonality>100_kg_m2"]))
diag.prec <- get.diag(df = DG, fac = factor(DG.meta[,"ECO_CHE_Bio12_prec.annual_kg>2500_kg_m2"]))

diag.eco.dry <- get.diag(df = DG, fac = factor(DG.meta[,"ECO_ecoregion_Dinerstein2017_Madagascar dry deciduous forests"]))
diag.eco.humid <- get.diag(df = DG, fac = factor(DG.meta[,"ECO_ecoregion_Dinerstein2017_Madagascar humid forests"]))
diag.eco.succulent <- get.diag(df = DG, fac = factor(DG.meta[,"ECO_ecoregion_Dinerstein2017_Madagascar succulent woodlands"]))

diag.eco.limestone <- get.diag(df = DG, fac = factor(DG.meta[,"ECO_surfaceLithologyAfrica_2 Limestone [Karst, Marls, Tsingy]"]))
diag.eco.sandstone <- get.diag(df = DG, fac = factor(DG.meta[,"ECO_surfaceLithologyAfrica_3 Sandstone [Non-Carbonate]"]))
diag.eco.basementrocks <- get.diag(df = DG, fac = factor(DG.meta[,"ECO_surfaceLithologyAfrica_4 Basement Rocks [Metasedimentary, Silicic, Metaigneous]"]))
diag.eco.alluvium <- get.diag(df = DG, fac = factor(DG.meta[,"ECO_surfaceLithologyAfrica_18 Alluvium [Various]"]))


res.diag <- rbind(diag.spec, diag.geo, diag.sg) # 909
best.diag <- subset(res.diag, TPR >= 0.9 & FPR <= 0.3) # 125 good ones

# write results
write.table(res.diag, file = paste0("Figures/DIAGNOSTIC_", nrow(DG), "_", ncol(DG),  "_", nrow(res.diag), ".txt"), quote = F, row.names = F, sep = "\t")
write.table(best.diag, file = paste0("Figures/DIAGNOSTIC_", nrow(DG), "_", ncol(DG),  "_", nrow(best.diag), ".txt"), quote = F, row.names = F, sep = "\t")

# best results
top.diag <- subset(res.diag, TPR >= 0.85 & FPR <= 0.15 & LEVEL <= 2) # 21 (52 at all levels, 11 at 1 level)


# VD and age
d.vd <- data.frame(dn[dn.meta$Species == "lemurica","VD",drop=F], dn.meta[dn.meta$Species == "lemurica","Species",drop=F], apply(data[rownames(dn[dn.meta$Species == "lemurica",]),c("Tree_height_m","Trunk_diameter_cm")], 2, as.numeric))
# d.vd <- data.frame(dn[,c("VD"),drop=F], dn.meta[,"Species",drop=F],apply(data[rownames(dn),c("Tree_height_m","Trunk_diameter_cm")], 2, as.numeric))


pvd1 <- ggplot(d.vd, aes(x = Trunk_diameter_cm, y = VD)) +
  geom_point() +
  geom_smooth() +
  labs(x = "Trunk diameter (cm)", y = expression(Vessel~density~(vessels/mm^{"2"}))) +
  facet_wrap(~Species) +
  theme_bw()

pvd2 <- ggplot(d.vd, aes(x = Tree_height_m, y = VD)) +
  geom_point() +
  geom_smooth() +
  labs(x = "Trunk height (m)", y = expression(Vessel~density~(vessels/mm^{"2"}))) +
  facet_wrap(~Species) +
  theme_bw()

ggpubr::ggarrange(pvd1, pvd2)

lm.vd <- lm(log(VD) ~ Tree_height_m + Trunk_diameter_cm, data = d.vd)
summary(lm.vd)

par(mfrow=c(2,2))
plot(lm.vd)
par(mfrow=c(1,1))

dsuc <- data[deco[deco$ECO_ecoregion_Dinerstein2017 %in% "Madagascar succulent woodlands","ID_Lab"],c("ID_Lab","Genus","Species")]
dsuc <- dsuc[dsuc$Genus %in% "Dalbergia",]
table(dsuc$Species) # brac, cherm, chloro, grev, lemu, perv, purp, rakotovaoi, tricho, tricolor, xero

# unique combinations
# WP_a & WP_c = 0 & WP_b = 1      Greveana   15/16        0

# boxplots
var <- "VL"
group <- "Species"
# group <- "ZT_Dalbergia_Group"
# group <- "Geogroup"
# thr <- 1 # low = "< thr" ; high = ">= thr"
thr <- 350


table(group = DG.meta[,group], variable = DG[,var] >= thr, useNA="a")

ggplot(cbind(DG, DG.meta[,!names(DG.meta) %in% names(DG)]),
       aes_string(x = group, y = var)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = thr)) +
  coord_flip() +
  theme_bw()


ggplot(cbind(DG, DG.meta[,!names(DG.meta) %in% names(DG)]),
       aes_string(x = group, y = var)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = thr)) +
  coord_flip() +
  facet_wrap(~Dataset) +
  theme_bw()


message("Do Table 2?")
scan()

###############
### Table 2 ###
###############

ids <- rownames(dn)
grp <- gsub(" s.l.", "", dn.meta$Species)
grp_id <- paste0("D_", grp, "_", ids)
ids.atlas <- gsub("Atlas_", "", rownames(da))

tab2 <- data.frame(Species = levels(dn.meta$Species_N))
# dnQ <- data.frame(apply(dn[,var.quant], 2, function(x) {rep(NA, length(x))}))
# rownames(dnQ) <- rownames(dn)
# dnQ$RW <- NA
dnQ <- dn[,0]

var.tvd <- c("TVD_RingPorousEW","TVD_RingPorousLW","TVD_DiffusePorous")
res <- SUM <- list()
for (i in c(var.quant, var.tvd)) {
  
  # df <- data.frame(readxl::read_excel("RawData/20220202_Dalbergia_quantitative_long.xlsx", sheet = i, .name_repair = "minimal"))
  df <- data.frame(readxl::read_excel(input, sheet = i, .name_repair = "minimal"))
  
  # grpQ <- sub("[.][0-9]+", "", names(df))
  # idsQ <- as.character(unname(df[1,]))
  grpQ <- sapply(strsplit(df$id, split = "_"), function(x) paste0(x[1], "_", x[2]))
  idsQ <- unlist(sapply(strsplit(df$id, split = "_"), "[", 3))
  grp_idQ <- paste0(grpQ, "_", idsQ)

  # check  
  stopifnot(!any(duplicated(idsQ))) # duplicated inds
  stopifnot(all(grp_id %in% grp_idQ)) # species_ind
  stopifnot(length(ids[!ids %in% idsQ]) == 0) # missing inds
  
  # extra inds
  extra <- idsQ[!idsQ %in% c(ids, ids.atlas)]
  
  # atlas inds
  atlas <- idsQ[idsQ %in% ids.atlas]

  # # read as 1-header table
  # dfQ <- data.frame(readxl::read_excel("RawData/20220202_Dalbergia_quantitative_long.xlsx", sheet = i, skip = 1, na = c("", " ", NA)))
  # stopifnot(all(sapply(names(dfQ), function(x) {class(dfQ[,x])}) %in% c("numeric","logical"))) # check numeric
  dfQ <- data.frame(t(df[,-1]))
  names(dfQ) <- idsQ
  
  # sort
  stopifnot(all(rownames(dn) %in% names(dfQ)))
  dfQ <- dfQ[,rownames(dn)]
  dfQ <- data.frame(apply(dfQ, 2, function(x) {c(x[!is.na(x)], x[is.na(x)])}))
  
  # compute individual median
  med <- apply(dfQ, 2, median, na.rm = TRUE)
  men <- apply(dfQ, 2, mean, na.rm = TRUE)
  n <- apply(dfQ, 2, function(x) {length(x[!is.na(x)])})
  
  # bind
  stopifnot(all.equal(rownames(dnQ), names(med)))
  dnQ[,paste0(i,"_median")] <- med
  dnQ[,paste0(i, "_mean")] <- men
  
  # compile data
  stopifnot(all.equal(names(med), rownames(dn)))
  if (!is.null(dn[names(med),i])) {
    d.comp <- data.frame(Ravo = dn[names(med),i], med = med, mean = men, n = n,
                         Species = dn.meta[names(med),"Species_N"],
                         grp = paste0("D_", gsub(" s.l.", "", dn.meta[names(med),"Species"]), "_", names(med)))
  } else {
    d.comp <- data.frame(Ravo = NA, med = med, mean = men, n = n,
                         Species = dn.meta[names(med),"Species_N"],
                         grp = paste0("D_", gsub(" s.l.", "", dn.meta[names(med),"Species"]), "_", names(med)))
  }
  # d.comp <- d.comp[!is.na(d.comp$Species),] # remove extra individuals
  stopifnot(nrow(d.comp) == nrow(dn)) # check if 93
  res[[i]] <- d.comp
  
  ## compare with Ravo data  
  # cor(d.comp$Ravo, d.comp$med, use = "pairwise")
  # plot(d.comp[,1:2]) ; abline(0,1)
  
  # # statistics
  # ggplot(data = d.comp, aes(x = med, y = Species)) +
  #   geom_boxplot() +
  #   # coord_flip() +
  #   theme_bw()
  # ggplot(data = d.comp, aes(x = mean, y = Species, fill = dn.meta$ZT_Dalbergia_Supergroup)) +
  #   geom_boxplot() +
  #   labs(x = paste0("mean ", i), fill = "Supergroup") +
  #   # coord_flip() +
  #   theme_bw()

  # # write transposed data to long format
  # dfQT <- data.frame(t(dfQ))
  # rownames(dfQT) <- d.comp$grp
  # write.table(dfQT, file = paste0("RawData/", i, ".csv"), sep = ";",
  #             quote = FALSE, row.names = TRUE, col.names = FALSE, na = "")
  
  ## compute species means, range of mean, SD of mean
  sp.mean <- tapply(d.comp$med, INDEX = d.comp$Species, FUN = mean, na.rm = T)
  sp.min <- tapply(d.comp$med, INDEX = d.comp$Species, FUN = min, na.rm = T)
  sp.max <- tapply(d.comp$med, INDEX = d.comp$Species, FUN = max, na.rm = T)
  sp.sd <- tapply(d.comp$med, INDEX = d.comp$Species, FUN = sd, na.rm = T)
  
  ## fill table
  if (mean(sp.mean,na.rm=T) > 100) r <- 0 else r <- 2
  cat("using", r, "digits\n")
  
  tab2[,i] <- paste0(round(sp.mean, r), " ± ", round(sp.sd, r), " (", round(sp.min, r), "–", round(sp.max, r), ")")
  tab2[,paste0("n_", i)] <- tapply(d.comp$med, INDEX = d.comp$Species, FUN = function(x) {sum(!is.na(x))})
  
  SUM[[i]] <- data.frame(Species = tab2$Species, Variable = i, n = tab2[,paste0("n_", i)],
                         min = sp.min, mean = sp.mean, sd = sp.sd, max = sp.max)
}

# write.table(tab2, file = "RawData/20220223_Table2_pretty.txt", quote = F, row.names = F, sep = "\t")
# write.table(dnQ, file = "RawData/20220223_Q.csv", quote = FALSE, row.names = TRUE, sep = ";")
# write.table(bind_rows(SUM), file = "RawData/20220223_Table2_numeric.txt", quote = F, row.names = F, sep = "\t")



##############################
### Test East vs. nonEast ####
##############################

# Mann-Whitney test (two-sample Wilcoxon test) on Ravo + Atlas dataset
dt <- data.frame(TVD_DiffusePorous = dc$TVD_DiffusePorous, TVD_RingPorousEW = dc$TVD_RingPorousEW, TVD_RingPorousLW = dc$TVD_RingPorousLW, VD = dc$VD, VL = dc$VL, dc.meta[,c("Species","ECO_ecoregion_Dinerstein2017_Madagascar humid forests")])
names(dt)[grep("humid",names(dt))] <- "humid"
dt[dt$Species == "orientalis" & dt$humid == 0,"humid"] <- 1
table(dt$Species, dt$humid)

dt2 <- dt[!dt$Species %in% c("lemurica","davidii"),]

# VD versus humid
aggregate(dt$VD, by = list(dt$humid), FUN = summary)
aggregate(dt$VD, by = list(dt$humid), FUN = sd)

wilcox.test(VD ~ humid, data = dt)
t.test(VD ~ humid, data = dt)

wilcox.test(VD ~ humid, data = dt2)
t.test(VD ~ humid, data = dt2)

aggregate(dt2$VD, by = list(dt2$humid), FUN = summary)
aggregate(dt2$VD, by = list(dt2$humid), FUN = sd)

# VL
wilcox.test(VL ~ humid, data = dt)


# TVD
summary(dc[dc$WP_c == 0,c("TVD_RingPorousEW","TVD_RingPorousLW")])
summary(dc[dc$WP_c == 1,c("TVD_DiffusePorous")])

aggregate(dc$TVD_RingPorousEW, by = list(diffuse = factor(dc$WP_c)), FUN = summary)
aggregate(dc$TVD_RingPorousLW, by = list(diffuse = factor(dc$WP_c)), FUN = summary)
aggregate(dc$TVD_DiffusePorous, by = list(diffuse = factor(dc$WP_c, levels = c("1","0"))), FUN = summary)


# TVD_Ring_PorousLW vs humid
aggregate(dt$TVD_DiffusePorous, by = list(dt$humid), FUN = summary)
aggregate(dt$TVD_DiffusePorous, by = list(dt$humid), FUN = sd, na.rm = T)

aggregate(dt$TVD_RingPorousEW, by = list(dt$humid), FUN = summary)
aggregate(dt$TVD_RingPorousEW, by = list(dt$humid), FUN = sd, na.rm = T)

aggregate(dt$TVD_RingPorousLW, by = list(dt$humid), FUN = summary)
aggregate(dt$TVD_RingPorousLW, by = list(dt$humid), FUN = sd, na.rm = T)

wilcox.test(TVD_RingPorousEW ~ humid, data = na.omit(dt[,c("TVD_RingPorousEW","humid")]))
wilcox.test(TVD_RingPorousLW ~ humid, data = na.omit(dt[,c("TVD_RingPorousLW","humid")]))
wilcox.test(TVD_DiffusePorous ~ humid, data = na.omit(dt[,c("TVD_DiffusePorous","humid")]))

# t.test(TVD_Ring_PorousLW ~ humid, data = na.omit(dt[,c("TVD_Ring_PorousLW","humid")]))



############################
### Supplementary Tables ###
############################

data$CODE <- data$ID_Lab
vars <- c("SpecimenID","CODE","Collection","Genus","Species",
          "CurrentDetermination","DeterminationBy",
          "MinimumElevation",
          # "MinimumDay","MinimumMonth","MinimumYear",
          "UpperName","LowerName")

### Appendix I: Collection details on Ravo 93
ds1 <- data[rownames(dn), vars]

### Appendix II: Raw data on Ravo 93
ds2 <- data.frame(merge(data[,c("CODE","Collection")], data.frame(CODE = rownames(dn), Genus = dn.meta$Genus, Species = dn.meta$Species, dn), by = "CODE", all.x = F, all.y = T, sort = F))
rownames(ds2) <- ds2$CODE
stopifnot(all(rownames(ds2) %in% rownames(ds1)))
ds2 <- data.frame(SpecimenID = ds1$SpecimenID, ds2[rownames(ds1),])
stopifnot(all.equal(ds1[,vars[1:4]], ds2[,vars[1:4]]))

### Appendix III: Identification Updates Atlas & Détienne 
ds3 <- data[grepl("Atlas|Détienne", data$data_wanatomy),c(vars[1:7], "PREVIOUS2", vars[8:length(vars)], "data_wanatomy")] # 85
names(ds3)[grep("data_wanatomy", names(ds3))] <- "Dataset"
names(ds3)[grep("PREVIOUS2", names(ds3))] <- "OldDetermination"
ds3 <- ds3[order(ds3$Dataset, sub("^cf._", "", ds3$Species)),]

### Appendix IV: 8 bioclimatic variables extracted from Karger: Ravo 93 + Atlas 8
# mean annual air temperature (Bio01)
# isothermality (i.e., the ratio of diurnal variation to annual variation in temperatures; Bio03)
# temperature seasonality (Bio04)
# annual range of air temperature (Bio07)
# annual precipitation amount (Bio12)
# precipitation amount of the wettest month (Bio13)
# precipitation amount of the driest month (Bio14)
# precipitation seasonality (Bio15)
var.clim <- c("ID_Lab",names(dc.meta)[grep("E_Bio01|E_Bio03|E_Bio04|E_Bio07|E_Bio12|E_Bio13|E_Bio14|E_Bio15",names(dc.meta))])
ds4 <- merge(data[,c("SpecimenID","ID_Lab","Collection","Genus","Species","data_wanatomy")], dc.meta[,var.clim], by = "ID_Lab", all.x = F, all.y = T, sort = F)
names(ds4) <- sub("(CHELSA_Bio\\d+)_([A-Za-z.]+)_(.*$)", "\\1_\\2 [\\3]", sub("kg_m2", "kg/m2", sub("degC", "°C", sub("^ECO_CHE_", "CHELSA_", names(ds4)))))
names(ds4)[grep("ID_Lab", names(ds4))] <- "CODE"
names(ds4)[grep("data_wanatomy", names(ds4))] <- "Dataset"
ds4 <- ds4[order(ds4$Dataset, ds4$Species, ds4$CODE), names(ds4)[c(2,1,3:ncol(ds4))]]

### Write tables
write.table(ds1, file = "Tables/TableS1.txt", row.names = F, quote = F, sep = "\t")
write.table(ds2, file = "Tables/TableS2.txt", row.names = F, quote = F, sep = "\t")
write.table(ds3, file = "Tables/TableS3.txt", row.names = F, quote = F, sep = "\t")
write.table(ds4, file = "Tables/TableS4.txt", row.names = F, quote = F, sep = "\t")

