## Helperfunctions
impute.na <- function(x, method = c("mean"), groups = NULL) {
  stopifnot(method %in% c("mean","median","zero","group.mean"))
  switch(method,
         mean = {x[is.na(x)] <- mean(x, na.rm = TRUE)},
         median = {x[is.na(x)] <- median(x, na.rm = TRUE)},
         zero = {x[is.na(x)] <- 0},
         group.mean = {
           stopifnot(!is.null(groups),
                     length(groups) == length(x))
           for (i in as.character(unique(groups))) {
             idx <- groups %in% i
             x[idx][is.na(x)[idx]] <- mean(x[idx], na.rm = TRUE)
           }
         }
  )
  return(x)
}
do.pairs <- function(df, method = "spearman", scale = TRUE, ellipses = TRUE) {
  # http://www.sthda.com/english/wiki/scatter-plot-matrices-r-base-graphs
  library(psych)
  pairs.panels(df, 
               method = "spearman",
               scale = scale,
               hist.col = "#00AFBB",
               density = TRUE,  # show density plots
               ellipses = ellipses # show correlation ellipses
  )
}
do.rpart <- function(df, fac, nmin = 5, control, plot = TRUE) {
  
  ## helperfunctions
  prune.rpart2 <- function(tree, cp) {
    ff <- tree$frame
    id <- as.integer(row.names(ff))
    toss <- id[ff$complexity <= cp & ff$var != "<leaf>"]
    if (length(toss) == 0L) 
      return(tree)
    newx <- snip.rpart(tree, toss)
    temp <- pmax(tree$cptable[, 1L], cp)
    keep <- match(unique(temp), temp)
    newx$cptable <- tree$cptable[keep, , drop = FALSE]
    newx$cptable[as.character(max(keep)), 1L] <- cp # fix
    # newx$variable.importance <- importance(newx) # hidden function
    newx
  }
  library(rpart)
  library(rpart.plot)
  library(adegenet)
  
  ## check input
  stopifnot(any(is.data.frame(df), is.matrix(df)),
            is.character(fac), fac %in% colnames(df), is.factor(df[,fac]),
            is.numeric(nmin), nmin > 0,
            is.list(control),
            is.logical(plot))
  
  ## prepare output
  factab <- table(df[,fac])
  if (nmin > sort(factab, decreasing = T)[2]) stop("nmin too large (only a single class remaining)!")
  levomit <- names(which(factab < nmin))
  if (length(levomit) > 0) df.omit <- df[-which(df[,fac] %in% levomit),] else df.omit <- df
  
  ## large rpart
  # set.seed(1234)
  d.train <- df.omit
  d.train[,fac] <- droplevels(d.train[,fac])
  parms <- list(prior = rep(1/nlevels(d.train[,fac]), nlevels(d.train[,fac])))
  rp.large <- rpart(as.formula(paste(fac, "~", ".")), data = d.train, control = control, parms = parms)
  
  ## complexity parameter pruning
  cp.tab <- rp.large$cptable
  cp.err <- cp.tab[,"xerror"]
  cp.sel <- cp.tab[min(which(cp.err <= cp.err[which.min(cp.err)] + 
                               cp.tab[,"xstd"][which.min(cp.err)])), "CP"]
  
  # rule: cp threshold for pruning at 1 SE above sd of minimum of the curve
  ethresh <- min(cp.tab[,"xerror"]) + cp.tab[which.min(cp.tab[,"xerror"]),"xstd"] 
  cp <- cp.tab[which(cp.tab[,"xerror"] < ethresh),"CP"][1]
  
  ## pruned rpart
  rp.final <- try(rpart::prune.rpart(rp.large, cp = cp), silent = TRUE)
  if (inherits(rp.final, "try-error")) {
    rp.final <- prune.rpart2(tree = rp.large, cp = cp)
  }

  ## plot pruned rpart model
  if (plot) {
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    split.fun <- function(x, labs, digits, varlen, faclen) {
      labs <- gsub("<", "\n<", gsub(">", "\n>", gsub(" ", "\n", gsub(":", "\n", labs))))
      for (i in 1:length(labs)) {
        
        # split labs[i] into multiple lines
        labs[i] <- paste(strwrap(labs[i], width = 10), collapse = "\n")
        rm(i)
      }
      labs
    }
    p <- capture.output(
      prp(rp.final, type = 3, extra = 0, # or use extra = 1
          box.col = transp(gg_color_hue(nlevels(df[,fac]))[rp.final$frame$yval], .6),  
          fallen.leaves = F, clip.right.labs = T, varlen = 0, space = 0, cex = 1.5, 
          gap = 8, ygap = 3, xflip = F, boxes.include.gap = T,
          under = T, branch.lty = 2, split.fun = split.fun)
    )
  }
  
  ## return results
  accuracy <- 1-sum(residuals(rp.final) == 1)/nrow(d.train)
  res <- list()
  res[["rp.input"]] <- d.train
  res[["rp.large"]] <- rp.large
  res[["cp.thresh"]] <- cp
  res[["rp.final"]] <- rp.final
  res[["accuracy"]] <- accuracy
  return(res)
}
layer_pcavar <- function(PCA, x = 1, y = 2, f = 1, quantile = 0.5, col = "tomato", cex = 1, length = 0.2) {
  
  library(ggplot2)
  library(ggrepel)
  
  # create scaled variable space (scaling) in discrimination space (space)
  fac <- diff(range(PCA[["xy"]][,c(x,y)]))/2.5
  m <- max(abs(PCA$rotation[,c(x,y)]))
  PCA$rotation.scaled <- apply(PCA$rotation[,c(x,y)], 2, function(x) x/(m/fac))
  
  # get arrow lengths
  l <- apply(PCA$rotation.scaled, 1, function(x) sqrt(x[1]^2 + x[2]^2))
  t <- quantile(l, quantile)
  a <- PCA$rotation.scaled[l >= t,]
  
  # only display top-quant arrows
  darr <- data.frame(x = a[,x]*f*1.05, y = a[,y]*f*1.05, row.names = rownames(a))
  arrows(x0 = 0, y0 = 0, x1 = a[,x]*f, y1 = a[,y]*f, 
         col = col, length = length, xpd = TRUE)
  text(a[,c(x,y)]*f*1.05, labels = rownames(a), 
       col = col, cex = cex, xpd = TRUE)
  
  # return repelled coordinates
  library(ggplot2)
  p <- ggplot(data.frame(PCA$xy), aes_string(x = colnames(PCA$xy)[1], y = colnames(PCA$xy)[2])) +
    geom_point() +
    # geom_segment(data = darr, x = 0, y = 0, aes(xend = x, yend = y),
    #              arrow = arrow(length = unit(0.15, "inches"))) +
    geom_text_repel(data = rbind(darr,darr), aes(x = x, y = y), label = c(rownames(darr), rep("", nrow(darr))), size = 2)
  
  p <- ggplot_build(p)$data[[2]][seq(nrow(darr)),1:2]
  rownames(p) <- rownames(darr)
  invisible(p)
}
layer_ldavar <- function(LDA, x = 1, y = 2, f = 1, quantile = 0.5, col = "tomato", cex = 1, length = 0.2) {
  
  library(ggplot2)
  library(ggrepel)
  
  # create scaled variable space (scaling) in discrimination space (space)
  LDA$mod$space <- LDA$x %*% LDA$mod$scaling
  fac <- diff(range(LDA$mod$space[,c(x,y)]))/2.5
  m <- max(abs(LDA$mod$scaling[,c(x,y)]))
  LDA$mod$scaling.scaled <- apply(LDA$mod$scaling[,c(x,y)], 2, function(x) x/(m/fac))
  
  # get arrow lengths
  l <- apply(LDA$mod$scaling.scaled, 1, function(x) sqrt(x[1]^2 + x[2]^2))
  t <- quantile(l, quantile)
  a <- LDA$mod$scaling.scaled[l >= t,]
  
  # only display top-quant arrows
  arrows(x0 = 0, y0 = 0, x1 = a[,1]*f, y1 = a[,2]*f, 
         col = col, length = length, xpd = TRUE)
  text(a[,c(1,2)]*f*1.05, labels = rownames(a), 
       col = col, cex = cex, xpd = TRUE)
}
cplot.phylo <- function(phy, df = NULL, type = "unrooted",
                        tip.label = TRUE, label.offset = 0.1, lab4ut = NULL,
                        edge.color = "black", palette.edge = NULL, edge.lwd = 1, edge.lty = 1,
                        tip.color = "black", palette.tip = NULL, tip.cex = NULL,
                        legend.edge = T, legend.tip = T,
                        legend.edge.pos = "topleft",legend.tip.pos = "bottomright",
                        legend.edge.cex = 1, legend.tip.cex = 1, ...) {
  
  stopifnot(inherits(phy, "phylo"))
  
  n <- length(phy$tip.label)
  
  # palettes
  funky <- function (n) {
    ramp <- grDevices::colorRampPalette(c("#A6CEE3","#1F78B4","#B2DF8A",
                                          "#33A02C","#FB9A99","#E31A1C",
                                          "#FDBF6F","#FF7F00","#CAB2D6",
                                          "#6A3D9A","#FFFF99","#B15928"))
    
    ramp(n)
  }
  eco <- function(n) {
    ramp <- grDevices::colorRampPalette(c("yellow","purple","green4",
                                          "pink","red","darkgreen","orange"))
    ramp(n)
  }
  viridis <- function(n) {
    ramp <- grDevices::colorRampPalette(c("#440154FF","#482173FF","#433E85FF",
                                          "#38598CFF","#2D708EFF","#25858EFF",
                                          "#1E9B8AFF","#2BB07FFF","#51C56AFF",
                                          "#85D54AFF","#C2DF23FF","#FDE725FF"))
    ramp(n)
  }
  plasma <- function(n) {
    ramp <- grDevices::colorRampPalette(c("#0D0887FF","#3E049CFF","#6300A7FF",
                                          "#8707A6FF","#A62098FF","#C03A83FF",
                                          "#D5546EFF","#E76F5AFF","#F58C46FF",
                                          "#FDAD32FF","#FCD225FF","#F0F921FF"))
    ramp(n)
  }
  spectre <- function(n) {
    ramp <- grDevices::colorRampPalette(c("#A6CEE3","#438EC0","#63A8A0","#98D277",
                                          "#E6E600","#B89B74","#F16667","#7A0403",
                                          "#B5B5B5","#FE982C","#3BA432","#C3AAD2",
                                          "#7D54A5","#00204D","#EAD27A", "#B15928"))
    ramp(n)
  }
  if (is.null(palette.edge)) palette.edge <- funky
  if (is.null(palette.tip)) palette.tip <- funky
  
  if (!is.null(df)) {
    stopifnot(inherits(df, c("data.frame","matrix")))
    
    # edge colour
    if (edge.color %in% colnames(df)) {
      if (!is.factor(df[,edge.color])) df[,edge.color] <- factor(df[,edge.color])
      col.ecol <- edge.color
      edge.color <- palette.edge(nlevels(df[,edge.color]))[df[phy$tip.label, edge.color]]
      edges <- sapply(phy$tip.label, function(x) which.edge(phy = phy, group = x))
      ecol <- rep(1, nrow(phy$edge))
      ecol[edges] <- edge.color
      edge.color <- ecol
    } else {
      col.ecol <- NULL
    }
    
    # tip colour
    if (tip.color %in% colnames(df)) {
      if (!is.factor(df[,tip.color])) df[,tip.color] <- factor(df[,tip.color])
      col.tcol <- tip.color
      tip.color <- palette.tip(nlevels(df[,tip.color]))[df[phy$tip.label, tip.color]]
    } else {
      col.tcol <- NULL
    }
  } else {
    col.ecol <- col.tcol <- NULL
  }
  
  # estimate a good tip.cex
  if (is.null(tip.cex)) {
    if (n > 40) {
      tip.cex <- 0.5
    } else {
      tip.cex <- 0.995 + (6.93E-03)*n + (-1.04E-03)*n^2 + (1.42E-05)*n^3
    }
  }
  
  # plot
  plot(phy, type = type, show.tip.label = tip.label, label.offset = label.offset,
       edge.color = edge.color, edge.width = edge.lwd, edge.lty = edge.lty, lab4ut = lab4ut,
       tip.color = tip.color, cex = tip.cex)
  
  # legend
  if (!is.null(col.ecol) & legend.edge) {
    legend(legend.edge.pos, bty = "n", pch = 22, cex = legend.tip.cex,
           legend = levels(df[,col.ecol]), title = col.ecol,
           pt.bg = palette.edge(nlevels(df[,col.ecol])), ...)
  }
  if (!is.null(col.tcol) & legend.tip & !identical(col.tcol, col.ecol)) {
    legend(legend.tip.pos, bty = "n", pch = 22, cex = legend.edge.cex,
           legend = levels(df[,col.tcol]), title = col.tcol,
           pt.bg = palette.edge(nlevels(df[,col.tcol])), ...)
  }
}
cplot <- function(X, x = 1, y = 2, flipx = 1, flipy = 1, zoom = 0,
                  variates = "X", add = NULL, drop = NULL,
                  .fillfac = NULL, .colfac = NULL, .hullfac = NULL, .spiderfac = NULL, .shapefac = NULL, .alphafac = NULL, .sizefac = NULL,
                  palette.fill = NULL, palette.col = NULL, palette.hull = NULL, palette.spider = NULL,
                  points = NULL, point.size = 2.5, point.shape = 21, point.alpha = 0.4, point.lwd = 1,
                  labels = FALSE, labels.cex = 2.5, labels.add = labels,
                  hull = !is.null(.hullfac), hull.concavity = 100, hull.alpha = 0.2, hull.expand = unit(.15,"cm"), hull.labels = FALSE, hull.labels.cex = 6, hull.legend = TRUE,
                  spider = !is.null(.spiderfac), spider.alpha = 0.8, spider.lwd = 0.5,
                  biplot = FALSE, loadings = NULL, quantile = NULL, f = 1, biplot.col = c("tomato","blue"), biplot.cex = 2.5, biplot.lwd = 0.25, biplot.alpha = 0.5,
                  plot = TRUE, legend.pos = "bottom", legend.size = 5, legend.spacing = unit(0, "cm"), legend.key.size = unit(0.25, "cm"), ...) {
  
  ## Usage
  # X       data.frame, matrix or projection object of class prcomp, lda, and various related classes (see list of projected)
  # x       numeric     variable or component to plot on x axis
  # y       numeric     variable or component to plot on y axis
  # flipx   numeric     constant to multiply with x (-1 for flip)
  # flipy   numeric     constant to multiply with y (-1 for flip)
  # zoom    numeric     zoom in (positive) or out (negative). Zero is no zoom.
  # variates character  character string denoting the scores component, if that is a list (e.g. for mixOmics)
  # .fillfac factor or numeric   Additional variables (numeric or factor) mapped to X
  # .colfac factor or numeric
  # .hullfac factor
  # .spiderfac factor
  # .shapefac factor
  # .alphafac numeric or factor
  # .sizefac numeric or factor
  # add     data.frame  new data to be projected to the same plot using the loadings matrix
  # drop    character  name of rows to drop from the plot
  # loadings character vector with name(s) of loading matrix or matrices
  # variates character with name of scores matrix
  
  # author: simon.crameri@usys.ethz.ch, Aug 2021
  library(ggplot2)
  library(ggrepel)
  library(ggforce)
  
  ## Check and re-organize input -----------------------------------------------
  
  # determine object class
  X.class <- sub("^pca$|^spca$|^ipca$|^sipca$|^mixo_pls$|^mixo_spls$|^mixo_plsda$|^mixo_splsda$", "mixOmics", class(X)[1])
  if ("prcomp" %in% class(X)) X.class <- "prcomp"
  if ("dudi" %in% class(X)) X.class <- "dudi"
  projected <- c("prcomp","princomp","Pca","dudi","Kpca","mixOmics","dapc","glPca","lda","LDA") # "Kplsr")
  
  # check input
  stopifnot(X.class %in% c("data.frame","matrix","umap","Kplsr","MDS","NMDS",projected))
  
  # factorize character vectors
  if (!is.null(.colfac)) if (is.character(.colfac)) .colfac <- as.factor(.colfac)
  if (!is.null(.fillfac)) if (is.character(.fillfac)) .fillfac <- as.factor(.fillfac)
  if (!is.null(.shapefac)) if (is.character(.shapefac)) .shapefac <- as.factor(.shapefac)
  if (!is.null(.alphafac)) if (is.character(.alphafac)) .alphafac <- as.factor(.alphafac)
  if (!is.null(.sizefac)) if (is.character(.sizefac)) .sizefac <- as.factor(.sizefac)
  if (!is.null(.hullfac)) if (is.character(.hullfac)) .hullfac <- as.factor(.hullfac)
  if (!is.null(.spiderfac)) if (is.character(.spiderfac)) .spiderfac <- as.factor(.spiderfac)
  
  # (re)set points argument
  if (is.null(points)) {
    points <- any(c(!is.null(.colfac), !is.null(.fillfac), !is.null(.shapefac), !is.null(.alphafac), !is.null(.sizefac)))
  }
  if (!points & all(c(is.null(.colfac), is.null(.fillfac), is.null(.shapefac), is.null(.alphafac), is.null(.sizefac), is.null(.hullfac), is.null(.spiderfac)))) {
    points <- TRUE 
  }
  
  # reset loadings argument
  if (is.null(loadings)) {
    loadings <- c("X","Y","var.load","loadings")
  }
  
  # reset palette argument
  funky <- function (n) {
    ramp <- grDevices::colorRampPalette(c("#A6CEE3","#1F78B4","#B2DF8A",
                                          "#33A02C","#FB9A99","#E31A1C",
                                          "#FDBF6F","#FF7F00","#CAB2D6",
                                          "#6A3D9A","#FFFF99","#B15928"))
    
    ramp(n)
  }
  eco <- function(n) {
    ramp <- grDevices::colorRampPalette(c("yellow","purple","green4",
                                          "pink","red","darkgreen","orange"))
    ramp(n)
  }
  viridis <- function(n) {
    ramp <- grDevices::colorRampPalette(c("#440154FF","#482173FF","#433E85FF",
                                          "#38598CFF","#2D708EFF","#25858EFF",
                                          "#1E9B8AFF","#2BB07FFF","#51C56AFF",
                                          "#85D54AFF","#C2DF23FF","#FDE725FF"))
    ramp(n)
  }
  plasma <- function(n) {
    ramp <- grDevices::colorRampPalette(c("#0D0887FF","#3E049CFF","#6300A7FF",
                                          "#8707A6FF","#A62098FF","#C03A83FF",
                                          "#D5546EFF","#E76F5AFF","#F58C46FF",
                                          "#FDAD32FF","#FCD225FF","#F0F921FF"))
    ramp(n)
  }
  spectre <- function(n) {
    ramp <- grDevices::colorRampPalette(c("#A6CEE3","#438EC0","#63A8A0","#98D277",
                                          "#E6E600","#B89B74","#F16667","#7A0403",
                                          "#B5B5B5","#FE982C","#3BA432","#C3AAD2",
                                          "#7D54A5","#00204D","#EAD27A", "#B15928"))
    ramp(n)
  }
  if (is.null(palette.spider)) palette.spider <- funky
  if (is.null(palette.hull)) palette.hull <- eco
  if (is.null(palette.col)) palette.col <- viridis
  if (is.null(palette.fill)) palette.fill <- plasma
  
  # specify internal parameters
  add_biplot <- biplot & X.class %in% projected
  add_to_plot <- !is.null(add) & X.class %in% projected
  drop_from_plot <- !is.null(drop)
  
  
  ## Exctract S (scores), L (loadings) and V (variance) comopnents -------------
  
  # scores matrix
  S <- switch(X.class,
              "data.frame" = data.frame(X),
              "matrix" = data.frame(X),
              "umap" = data.frame(X$layout),
              "prcomp" = data.frame(X$x),
              "princomp" = data.frame(X$scores),
              "Pca" = data.frame(X$T),
              "Kpca" = data.frame(X$T),
              "dudi" = data.frame(X$li),
              "mixOmics" = data.frame(X$variates[[variates]]), # $x is same as $variates$X
              "dapc" = data.frame(X$ind.coord),
              "glPca" = data.frame(X$scores),
              "lda" = data.frame(predict(X)$x),
              "LDA" = data.frame(X$mod.pred$x),
              "Kplsr" = data.frame(X$T),
              "MDS" = data.frame(X$x),
              "NMDS" = data.frame(X$points))
  
  # loadings matrix
  L <- switch(X.class,
              # "data.frame" = NULL,
              # "matrix" = NULL,
              # "umap" = NULL,
              "prcomp" = X$rotation,
              "princomp" = X$loadings,
              "Pca" = X$P,
              "Kpca" = X$P,
              "dudi" = {
                L <- X$c1
                # colnames(L) <- colnames(S)
                as.matrix(L)
              },
              "mixOmics" = X$loadings[loadings], # $rotation is same as $loadings$X
              "dapc" = {
                L.pc <- !is.null(X$loadings)
                L.v <- !is.null(X$var.load)
                l.pc <- "loadings" %in% loadings
                l.v <- "var.load" %in% loadings
                
                # PCA loadings
                if (L.pc & l.pc & !l.v) L <- X$loadings
                
                # var loadings
                if (L.v & l.v & !l.pc) L <- X$var.load
                
                # PCA and var loadings
                if (L.v & L.pc & l.pc & l.v) L <- X[c("loadings","var.load")]
                
                L
                # if (!is.null(X$var.load)) X$var.load else X$loadings # X$loadings gives PCA loadings
              },
              "glPca" = {
                if (!is.null(X$loadings)) L <- X$loadings else L <- NULL
              },
              "lda" = X$scaling,
              "LDA" = X$mod$scaling,
              "Kplsr" = NULL,
              "MDS" = NULL,
              "NMDS" = NULL)
  
  # convert L matrix to list for compatibility with X and Y loadings (mixOmics)
  if (!is.list(L)) {
    L <- list(L)
    loadings <- 1
  } else {
    loadings <- loadings[loadings %in% names(L)]
  }
  
  # percent explained variation
  v <- switch(X.class,
              "data.frame" = NULL,
              "matrix" = NULL,
              "umap" = NULL,
              "prcomp" = X$sdev^2 / sum(X$sdev^2),
              "princomp" = X$sdev^2 / sum(X$sdev^2),
              "Pca" = X$eig / sum(X$eig),
              "dudi" = X$eig / sum(X$eig),
              "Kpca" = X$eig / sum(X$eig),
              "mixOmics" = X$prop_expl_var[[variates]],
              "dapc" = X$eig / sum(X$eig),
              "glPca" = X$eig / sum(X$eig),
              "lda" = X$svd^2 / sum(X$svd^2),
              "LDA" = X$mod$svd^2 / sum(X$mod$svd^2),
              "Kplsr" = NULL,
              "MDS" = NULL,
              "NMDS" = NULL)
  
  # component name
  n <- switch(X.class,
              "data.frame" = colnames(X),
              "matrix" = colnames(X),
              "umap" = "Comp.",
              "prcomp" = "PC",
              "princomp" = "PC",
              "Pca" = "PC",
              "dudi" = sub("[0-9 ]+", "", colnames(S)),
              "Kpca" = "KPC",
              "mixOmics" = unique(sub("[0-9 ]+", "", colnames(S))),
              "dapc" = "LD",
              "glPca" = "PC",
              "lda" = "LD",
              "LDA" = "LD",
              "Kplsr" = "LV",
              "MDS" = "MDS",
              "NMDS" = "MDS")
  
  compname <- function(X, x, class = X.class, v = NULL) {
    if (class %in% c("data.frame","matrix")) {
      n[x]
    } else if (class %in% c("umap","Kplsr","MDS","NMDS")) {
      paste(n, x)
    } else if (class %in% projected) {
      paste0(n, " ", x, " (", round(100*v[x], 2), "%)")
    }
  }
  
  ## Add or drop points from data ----------------------------------------------
  
  # add
  if (add_to_plot) {
    stopifnot(inherits(add, c("matrix","data.frame")),
              all(rownames(L[[loadings]]) %in% colnames(add)))
    
    Sadd <- as.matrix(add[,rownames(L[[loadings]])]) %*% L[[loadings]]
    Sadd[,x] <- Sadd[,x]*flipx
    Sadd[,y] <- Sadd[,y]*flipy
    Sadd <- data.frame(Sadd)
  }
  
  # drop
  if (drop_from_plot) {
    stopifnot(inherits(drop, c("matrix","data.frame","character","factor")))
    
    if (class(drop) %in% c("character","factor")) {
      if (!all(as.character(drop) %in% rownames(S))) {
        warning("Undefined rows selected in <drop>")
      }
      keep <- !rownames(S) %in% as.character(drop)
      S <- S[keep,]
    } else if (class(drop) %in% c("matrix","data.frame")) {
      if (!all(rownames(drop) %in% rownames(S))) {
        warning("Not all rownames(drop) in rownames(X)!")
      }
      keep <- !rownames(S) %in% rownames(drop)
      S <- S[keep,]
    } else if (class(drop) %in% c("integer")) {
      if (!all(drop) %in% seq_along(nrow(S))) {
        warning("Undefined rows selected in <drop>")
      }
      keep <- !seq_along(nrow(X)) %in% drop
      S <- S[keep,]
    }
    .hullfac <- .hullfac[keep]
    .spiderfac <- .spiderfac[keep]
    .colfac <- .colfac[keep]
    .fillfac <- .fillfac[keep]
    .shapefac <- .shapefac[keep]
    .alphafac <- .alphafac[keep]
    .sizefac <- .sizefac[keep]
  }
  
  ## Plot data -----------------------------------------------------------------
  
  # flip
  if (ncol(S) < 2) x <- y <- 1
  if (is.character(x) & x %in% colnames(S)) x <- which(colnames(S) == x)
  if (is.character(y) & y %in% colnames(S)) y <- which(colnames(S) == y)
  S[,x] <- S[,x]*flipx
  S[,y] <- S[,y]*flipy
  
  # visualize
  xlim <- range(S[,x])
  ylim <- range(S[,y])
  if (!is.null(.hullfac)) {
    fillname <- substitute(.hullfac)
  }  else if (!is.null(.fillfac)) {
    fillname <- substitute(.fillfac)
  } else {
    fillname = ""
  }
  if (!is.null(.spiderfac)) {
    colorname <- substitute(.spiderfac)
  }  else if (!is.null(.colfac)) {
    colorname <- substitute(.colfac)
  } else {
    colorname = ""
  }
  p <- ggplot(data = S, aes_string(x = colnames(S)[x], y = colnames(S)[y])) +
    labs(x = compname(X = X, x = x, v = v),
         y = compname(X = X, x = y, v = v),
         fill = fillname,
         color = "", alpha = "", size = "")
  
  if (hull) {
    
    col.hull <- palette.hull(nlevels(.hullfac))[which(levels(.hullfac) %in% unique(.hullfac))]
    notna <- rowSums(is.na(S[,c(x,y)])) == 0
    
    p <- p +
      geom_mark_hull(data = na.omit(S[,c(x,y)]),
                     inherit.aes = F, 
                     aes_string(x = colnames(S)[x], y = colnames(S)[y],
                                fill = .hullfac[notna], label = if (hull.labels) .hullfac[notna] else NULL),
                     concavity = hull.concavity, 
                     radius = unit(.15,"cm"),
                     expand = hull.expand,
                     alpha = hull.alpha,
                     size = .1,
                     label.fontsize = hull.labels.cex,
                     con.size = 0.25,
                     con.type = "elbow",
                     con.cap = unit(0.5, "mm"),
                     show.legend = hull.legend,
                     na.rm = TRUE)
    
    if (hull.legend) {
      p <- p +
        scale_fill_manual(values = col.hull) +
        guides(fill = guide_legend(title = substitute(.hullfac), title.theme = element_text(size = legend.size),
                                   # label.position = "bottom",
                                   override.aes = list(stroke = 0, size = 0, linetype = 0), # does not work!?
                                   title.position = "top")) # this is ignored?!
    }
    
  }
  
  if (spider) {
    if (is.factor(.spiderfac)) {
      centroids <- data.frame(apply(X = S[,c(x,y)], MARGIN = 2, FUN = function(x) {tapply(x, INDEX = droplevels(.spiderfac), FUN = mean, na.rm = TRUE)}))
      colnames(centroids) <- paste0(colnames(centroids), ".c")
      centroids$.spiderfac <- levels(droplevels(.spiderfac))
      Sc <- merge(data.frame(S[,c(x,y)], .spiderfac, id = seq(nrow(S))), centroids, by = ".spiderfac", all = T, sort = F)
      Sc <- Sc[order(Sc$id),]
    } else {
      Sc <- data.frame(rep(mean(S[,x],na.rm=T), nrow(S)), rep(mean(S[,y], na.rm=T), nrow(S)))
      colnames(Sc) <- paste0(colnames(S[,c(x,y)]), ".c")
      Sc <- cbind(S[,c(x,y)], Sc)
    }
    
    col.spider <- palette.spider(nlevels(.spiderfac))[which(levels(.spiderfac) %in% unique(.spiderfac))]
    
    p <- p + geom_segment(data = Sc, inherit.aes = FALSE,
                          aes_string(x = paste0(colnames(S)[x], ".c"), 
                                     xend = colnames(S)[x], 
                                     y = paste0(colnames(S)[y], ".c"),
                                     yend = colnames(S)[y],
                                     color = .spiderfac),
                          size = spider.lwd, alpha = spider.alpha,
                          show.legend = TRUE) +
      scale_color_manual(values = col.spider)
  }
  
  if (points) {
    # cleg = "legend"
    if (!is.null(.shapefac)) point.shape <- NULL
    
    if (!is.null(.colfac) & !is.null(.spiderfac) | !is.null(.fillfac) & !is.null(.hullfac)) {
      need.pckg <- "ggnewscale"
      if (any(!need.pckg %in% installed.packages())) {
        for (i in need.pckg[!need.pckg %in% installed.packages()]) {
          msg <- paste0("The <", i, "> package is needed. Do you wish to install it? [y/n]")
          q <- readline(prompt = message(msg))
          if (grepl("^y|^Y", q)) install.packages(i)
        }
      }
      if (!all(need.pckg %in% installed.packages())) {
        warning("Multiple color aestetics needed. Please install <ggnewscale>")
      }
    }
    if (!is.null(.colfac) & !is.null(.spiderfac)) {
      p <- p + 
        ggnewscale::new_scale_color() + labs(color = "")
    }
    if (!is.null(.fillfac) & !is.null(.hullfac)) {
      p <- p + 
        ggnewscale::new_scale_fill() + labs(fill = "")
    }
    
    # set dots (point.size only if .sizefac = NULL, etc. for alpha and shape)
    # dots <- list(...)
    
    if (!is.null(.fillfac) & is.null(.shapefac)) {
      p <- p +
        geom_point(inherit.aes = FALSE,
                   aes_string(x = colnames(S)[x], y = colnames(S)[y], 
                              color = .colfac, fill = .fillfac, size = .sizefac,
                              shape = .shapefac, alpha = .alphafac),
                   shape = point.shape, ...)
    } else {
      p <- p +
        geom_point(inherit.aes = FALSE,
                   aes_string(x = colnames(S)[x], y = colnames(S)[y], 
                              color = .colfac, fill = .fillfac, size = .sizefac,
                              shape = .shapefac, alpha = .alphafac),
                   ...)
    }
    
    # set scales
    if (!is.null(.shapefac)) {
      p <- p +
        scale_shape_manual(values = 21:25)
    }
    if (!is.null(.colfac) & is.factor(.colfac)) {
      col.color <- palette.col(nlevels(.colfac))[which(levels(.colfac) %in% unique(.colfac))]
      p <- p +
        scale_color_manual(values = col.color) # , guide = cleg
    }
    if (!is.null(.colfac) & is.numeric(.colfac)) {
      N <- 3
      col.color <- palette.col(N)
      p <- p +
        scale_color_gradient2(low = col.color[1], mid = col.color[(N+1)/2], high = col.color[N],
                              midpoint = mean(.colfac, na.rm = TRUE))
    }
    if (!is.null(.fillfac) & is.factor(.fillfac)) {
      col.fill <- palette.fill(nlevels(.fillfac))[which(levels(.fillfac) %in% unique(.fillfac))]
      p <- p +
        scale_fill_manual(values = col.fill)
    }
    if (!is.null(.fillfac) & is.numeric(.fillfac)) {
      N <- 3
      col.fill <- palette.fill(N)
      p <- p +
        scale_fill_gradient2(low = col.fill[1], mid = col.fill[(N+1)/2], high = col.fill[N],
                             midpoint = mean(.fillfac, na.rm = TRUE))
    }
  }
  
  if (add_biplot) {
    add.biplot <- function(L, quantile = NULL, f = 1) {
      # set quantile
      if (is.null(quantile)) {
        if (nrow(L) > 20) quantile <- 0.9 else quantile <- 0
      }
      
      # scale loadings to scores
      sfac <- diff(range(S[,c(x,y)]))/2.5
      m <- max(abs(L[,c(x,y)]))
      L.scaled <- apply(L[,c(x,y)], 2, function(x) x/(m/sfac))
      
      # get arrow lengths
      l <- apply(L.scaled, 1, function(x) sqrt(x[1]^2 + x[2]^2))
      t <- quantile(l, quantile)
      a <- L.scaled[l >= t,,drop=F]
      
      # subset top-quantile arrows
      darr <- data.frame(x = a[,1]*f, y = a[,2]*f, row.names = rownames(a))
    }
    
    # check quantile, f, biplot.cex, biplot.col
    if (length(loadings) > 1) {
      if (!is.null(quantile) & length(quantile) == 1) quantile <- rep(quantile, 2)
      if (length(f) == 1) f <- rep(f, 2)
      if (length(biplot.cex) == 1) biplot.cex <- rep(biplot.cex, 2)
      if (length(biplot.col) == 1) biplot.col <- rep(biplot.col, 2)
      if (length(biplot.lwd) == 1) biplot.lwd <- rep(biplot.lwd, 2)
      if (length(biplot.alpha) == 1) biplot.alpha <- rep(biplot.alpha, 2)
    }
    
    for (i in loadings) {
      
      id <- which(loadings == i)
      
      # create arrow segments
      darr <- add.biplot(L = L[[i]], quantile = quantile[id], f = f[id])
      darr$label <- rownames(darr)
      assign(paste0("darr.", i), darr)
      
      # plot arrows
      p <- p +
        geom_segment(data = get(paste0("darr.", i)), inherit.aes = FALSE,
                     aes(x = 0, y = 0, xend = x, yend = y),
                     arrow = arrow(angle = 30, length = unit(.25,"cm")),
                     size = biplot.lwd[id], alpha = biplot.alpha[id],
                     color = biplot.col[id]) +
        geom_text_repel(data = get(paste0("darr.", i)), inherit.aes = FALSE,
                        aes(x = x, y = y, label = label),
                        color =  biplot.col[id],
                        size = biplot.cex[id])
      xlim <- range(c(xlim, darr$x))
      ylim <- range(c(ylim, darr$y))
    }
  }
  
  if (add_to_plot) {
    p <- p +
      geom_point(data = Sadd, inherit.aes = FALSE,
                 aes_string(x = paste0(n, x), y = paste0(n, y)),
                 size = point.size*1.5, alpha = point.alpha*1.5,
                 fill = "tomato", colour = "black", shape = point.shape)
    xlim <- range(c(xlim, Sadd[,x]))
    ylim <- range(c(ylim, Sadd[,y]))
  }
  
  if (labels) {
    p <- p +
      geom_text_repel(inherit.aes = FALSE,
                      aes_string(x = paste0(n, x), y = paste0(n, y)),
                      label = rownames(S), size = labels.cex, segment.size = 0.2, ...)
  }
  
  if (add_to_plot & labels.add) {
    p <- p +
      geom_text_repel(data = Sadd, inherit.aes = FALSE,
                      aes_string(x = paste0(n, x), y = paste0(n, y)),
                      label = rownames(Sadd), size = labels.cex*1.2, segment.size = 0.2, ...)
  }
  
  # add zoom and theme
  margin <- c(-1,1)*zoom*-0.1
  xlim <- xlim + margin*abs(diff(xlim))
  ylim <- ylim + margin*abs(diff(ylim))
  
  p <- p +
    lims(x = xlim, y = ylim) +
    guides(color = guide_legend(title = substitute(.colfac), title.theme = element_text(size = legend.size),
                                label.position = "bottom",
                                override.aes = list(linetype = 0),
                                title.position = "top"),
           # fill = guide_legend(title = substitute(.fillfac), title.theme = element_text(size = legend.size),
           #                     label.position = "bottom",
           #                     override.aes = list(linetype = 0, shape = 21),
           #                     title.position = "top"),
           alpha = guide_legend(title = substitute(.alphafac), title.theme = element_text(size = legend.size),
                                label.position = "bottom",
                                override.aes = list(linetype = 0),
                                title.position = "top"),
           size = guide_legend(title = substitute(.sizefac), title.theme = element_text(size = legend.size),
                               label.position = "bottom",
                               override.aes = list(linetype = 0),
                               title.position = "top"),
           shape = guide_legend(title = substitute(.shapefac), title.theme = element_text(size = legend.size),
                                label.position = "bottom",
                                override.aes = list(linetype = 0),
                                title.position = "top")) +
    theme_bw() +
    theme(legend.position = legend.pos,
          legend.text = element_text(size = legend.size),
          legend.spacing = legend.spacing,
          legend.key.height = legend.key.size,
          legend.key.width = legend.key.size)
  
  # plot
  if (plot) {try(print(p), silent = TRUE)}
  
  # return
  invisible(p)
  
}
annotate.n <- function(x, valid = TRUE) {
  char <- as.character(x)
  tab <- table(char)
  if (!valid) {
    new <- sapply(unique(char), function(x) {char[char == x] <- paste0(x, " (n = ", tab[x], "]")})
  } else {
    new <- sapply(unique(char), function(x) {char[char == x] <- paste0(x, "_n", tab[x])})
  }
  newchar <- sapply(char, function(x) {new[x]})
  as.character(newchar)
}

###########################
### Tidymodels wrappers ###
###########################

tidytune <- function(Xtrain, Ytrain, Xtest, Ytest, workflow, tuneGrid, folds,
                     trControl, metrics, metric.best = "accuracy",
                     seed.train = NULL, seed.final = NULL,
                     plot.tune = TRUE, plot.roc = FALSE, ...) {
  
  # author: simon.crameri@usys.ethz.ch, Sep 2021
  
  t1 <- Sys.time()
  library(tidymodels)
  
  ## Check input
  stopifnot(inherits(workflow, "workflow"),
            inherits(tuneGrid, c("data.frame","matrix")),
            inherits(folds, "rset"),
            inherits(trControl, "control_resamples"),
            inherits(metrics, "metric_set"),
            metric.best %in% unlist(tibble(metrics = deparse(metrics)) %>%
                                      filter(grepl("UseMethod", metrics)) %>%
                                      mutate(metrics = sub('"\\)', '', 
                                                           sub('[ ]+UseMethod\\("','', 
                                                               metrics)))),
            is.logical(plot.tune), is.logical(plot.roc))
  
  ## Set parameters
  rec <- workflow$pre$actions$recipe$recipe
  mod <- workflow$fit$actions$model$spec
  
  y <- filter(rec$var_info, role == "outcome")$variable
  fform <- class(mod)[1]
  mode <- mod$mode
  engine <- mod$engine
  transl <- translate(mod)
  pkg_fun <- paste(transl$method$fit$func, collapse = "_")
  
  met_meta <- c(".metric",".estimator","mean","n","std_err",".config")
  
  ## Model fitting
  totune <- nrow(tuneGrid) > 1 | any(sapply(mod$args, as_label) %in% "tune()")
  if (totune) {
    
    ## Tune hyperparameters using tuneGrid, folds, metrics, trControl
    if (!is.null(seed.train)) {
      .Random.seed <- seed.train
    } else {
      seed.train <- .Random.seed
    }
    fit_train <- try(
      workflow %>% 
        tune_grid(resamples = folds, grid = tuneGrid,
                  metrics = metrics, control = trControl),
      silent = TRUE)
    attr(fit_train, "seed") <- seed.train # use .Random.seed <- attr(fit_train, "seed") to reproduce
    
    if (inherits(fit_train, "try-error")) {
      err <- fit_train[1]
      need.pkg <- grepl("Some package installs are required", err)
      
      # install needed package(s)
      if (need.pkg) {
        need.pkg <- gsub("'\\n|^'", '', sub('[A-Za-z :]+', '', err))
        for (i in need.pkg) {
          message("installing package <", i, ">")
          install.packages(i)
        }
      }
      
      # tune hyperparameters
      if (!is.null(seed.train)) {
        .Random.seed <- seed.train
      } else {
        seed.train <- .Random.seed
      }
      fit_train <- workflow %>% 
        tune_grid(resamples = folds, grid = tuneGrid,
                  metrics = metrics, control = trControl)
      attr(fit_train, "seed") <- seed.train # use .Random.seed <- attr(fit_train, "seed") to reproduce
    }
    
    ## Butcher data in splits component
    for (i in seq(length(fit_train$splits))) {
      fit_train$splits[[i]]$data <- tibble(.rows = nrow(fit_train$splits[[i]]$data))
    }
    
    ## Get best hyperparameter combination (consider select_by*)
    hpar <- colnames(tuneGrid)
    met_train <- collect_metrics(fit_train) %>% # summarizes fit_train$.metric elements
      dplyr::select(any_of(c(hpar, met_meta)))
    
    bestTune <- try(fit_train %>% select_best(metric.best) %>%
                      dplyr::select(any_of(c(hpar, met_meta))),
                    silent = TRUE)
    if (inherits(bestTune, "try-error")) {
      on.exit(return(fit_train))
      stop() 
    }
    
    ## Define final workflow
    wflow_final <- workflow %>% finalize_workflow(bestTune)
    
  } else {
    
    message("Nothing to tune, fitting full training data to workflow")
    
    ## Define final workflow
    fit_train <- bestTune <- met_train <- as_tibble(NA)
    wflow_final <- workflow
    
  }
  
  ## Fit full model with best hyperparameters
  s <- list(data = tibble(rbind(cbind(Xtrain, .Y = Ytrain),
                                cbind(Xtest, .Y = Ytest))) %>%
              rename(!!y := .Y),
            in_id = 1:nrow(Xtrain),
            out_id = (nrow(Xtrain)+1):(nrow(Xtrain)+nrow(Xtest)),
            id = tibble(id = "Resample1"))
  s <- structure(s, class = c("mc_split","rsplit"))
  
  if (!is.null(seed.final)) {
    .Random.seed <- seed.final
  } else {
    seed.final <- .Random.seed
  }
  fit_final <- last_fit(wflow_final, split = s, metrics = metrics)
  attr(fit_final$.workflow[[1]], "seed") <- seed.final # use .Random.seed <- attr(mod, "seed") to reproduce
  
  ## Predictions
  # training dataset
  if (mode == "classification") {
    pred.train <- tibble(
      predict(fit_final$.workflow[[1]], new_data = training(s), type = "class"),
      predict(fit_final$.workflow[[1]], new_data = training(s), type = "prob"),
      !!y := Ytrain)
  } else {
    pred.train <- predict(fit_final$.workflow[[1]], new_data = training(s))
  }
  
  # test dataset
  if (mode == "classification") {
    pred.test <- tibble(
      predict(fit_final$.workflow[[1]], new_data = testing(s), type = "class"),
      predict(fit_final$.workflow[[1]], new_data = testing(s), type = "prob"),
      !!y := Ytest)
  } else {
    # returns default prediction type
    pred.test <- fit_final$.predictions[[1]] %>% dplyr::select(-c(".row",".config"))
  }
  
  ## Performance
  met_final <- fit_final %>% collect_metrics() # also in fit_final$.metrics[[1]]
  
  # plot ROC
  predicted <- levels(droplevels(pred.test[[y]]))
  roc <- try(pred.test %>%
               roc_curve(truth = !!y,
                         estimate = !!(paste(paste0(".pred_", predicted)))),
             silent = TRUE)
  
  ## Butcher data in splits component
  for (i in seq(length(fit_final$splits))) {
    fit_final$splits[[i]]$data <- tibble(.rows = nrow(fit_final$splits[[i]]$data))
  }
  # for (i in seq(length(folds$splits))) {
  #   folds$splits[[i]]$data <- tibble(.rows = nrow(folds$splits[[i]]$data))
  # }
  
  ## List of arguments
  args <- list(Xtrain = tibble(Xtrain), Ytrain = Ytrain, Xtest = tibble(Xtest), Ytest = Ytest,
               # workflow = workflow, # trained workflow stored in res$mod
               # folds = folds, # fold$splits stored in res$train$splits (data removed)
               tuneGrid = tibble(tuneGrid), trControl = trControl,
               metrics = metrics, metric.best = metric.best,
               plot.tune = plot.tune, plot.roc = plot.roc, dots = list(...))
  
  ## Compile results
  nTT <- sapply(list(train = s$in_id, test = s$out_id,
                     total = c(s$in_id, s$out_id)), length)
  t2 <- Sys.time()
  res <- list(call = match.call(), mode = mode, fform = fform, engine = engine, func = pkg_fun,
              train = fit_train, nX = ncol(fit_final$.workflow[[1]]$pre$mold$predictors),
              nSplit = tibble(data.frame(set = names(nTT), n = nTT)),
              nTune = nrow(tuneGrid), bestTune = bestTune, 
              mod = fit_final$.workflow[[1]], notes = fit_final$.notes, 
              perf.train = met_train, perf.test = met_final,
              pred.train = pred.train, pred.test = pred.test,
              args = args, time = list(t1 = t1, t2 = t2, elapsed = t2 - t1,
                                       parallel = try(c(foreach::getDoParRegistered()), silent = T),
                                       name = try(c(foreach::getDoParName()), silent = T),
                                       ncpu = try(foreach::getDoParWorkers(), silent = T)))
  
  res <- structure(res, class = "tidytune")
  
  ## Plot results
  if (totune & plot.tune & !inherits(try(plot(res), silent = TRUE), "try-error")) {print(plot(res))}
  if (totune & plot.roc & !inherits(roc, "try-error")) {print(autoplot(roc))}
  
  ## Return
  invisible(res)
}

# combine training of two tidytune objects and re-fit the best model
c.tidytune <- function(res1, res2, plot.tune = TRUE) {
  
  stopifnot(inherits(res2, "tidytune"),
            inherits(res1, "tidytune"))
  
  t1 <- Sys.time()
  library(tidymodels)
  
  # check input
  stopifnot(inherits(res1, "tidytune"),
            inherits(res2, "tidytune"),
            identical(res1$mode, res2$mode),
            res1$fform == res2$fform,
            res1$engine == res2$engine,
            res1$func == res2$func)
  
  # set par
  metric.best <- unique(c(res1$args$metric.best, res2$args$metric.best)) # only first considered
  mfun <- ifelse(attr(eval(parse(text = metric.best[1])),"direction") == "maximize", max, min)
  
  met_meta <- c(".metric",".estimator","mean","n","std_err",".config")
  
  # bind tuning metrics
  hpar <- unique(c(names(res1$args$tuneGrid),
                   names(res2$args$tuneGrid),
                   ".set"))
  
  met_train <- bind_rows(mutate(res1$perf.train, .set = "1"),
                         mutate(res2$perf.train, .set = "2")) %>%
    dplyr::select(any_of(c(hpar, met_meta)))
  
  # select best
  bestTune <- met_train %>% 
    filter(.metric == metric.best) %>%
    filter(mean == mfun(mean)) %>% "["(1,)
  
  bestSet <- bestTune[[".set"]]
  otherSet <- c("1","2")[!c("1","2") %in% bestSet]
  best <- get(paste0("res", bestSet))
  
  # create tidytune structure based on combined results
  fit_train <- list(splits = NA,
                    id = NA,
                    id2 = NA,
                    .metrics = met_train,
                    .notes = NA)
  tuneGrid <- met_train %>% filter(.metric == metric.best) %>% 
    dplyr::select(-any_of(met_meta))
  nTune <- nrow(tuneGrid)
  args <- best$args
  args$tuneGrid <- tuneGrid
  
  t2 <- Sys.time()
  res <- list(call = match.call(), mode = best$mode, fform = best$fform, engine = best$engine, func = best$func,
              train = fit_train, nX = ncol(best$mod$pre$mold$predictors),
              nSplit = best$nSplit, nTune = nTune, bestTune = bestTune, 
              mod = best$mod, notes = best$notes, 
              perf.train = met_train, perf.test = best$perf.test,
              pred.train = best$pred.train, pred.test = best$pred.test,
              args = args, time = list(t1 = t1, t2 = t2, elapsed = t2 - t1,
                                       parallel = try(c(foreach::getDoParRegistered()), silent = T),
                                       name = try(c(foreach::getDoParName()), silent = T),
                                       ncpu = try(foreach::getDoParWorkers(), silent = T)))
  
  res <- structure(res, class = "tidytune")
  
  if (plot.tune) try(plot(res), silent = TRUE)
  invisible(res)
}

# return a summary string for file names
prefix.tidytune <- function(tidytune) {
  
  stopifnot(inherits(tidytune, "tidytune"))
  
  paste(tidytune$fform, tidytune$func, nrow(tidytune$args$Xtrain),
        table(tidytune$mod$pre$actions$recipe$recipe$var_info$role)["predictor"],
        nlevels(tidytune$args$Ytrain), nrow(tidytune$args$tuneGrid), sep = "_")
}

# tidytune plot methods: plot.tidytune
plot.tidytune <- function(tidytune, order = NULL,
                          x.log10 = TRUE, y.log10 = FALSE, mode = "alpha",
                          point.size = 1, point.alpha = 0.6, line.alpha = 0.6) {
  
  stopifnot(inherits(tidytune, "tidytune"))
  library(ggplot2)
  
  ## Get components
  tuneGrid <- tidytune$args$tuneGrid
  met_train <- tidytune$perf.train
  bestTune <- tidytune$bestTune
  pkg_fun <- tidytune$func
  
  ## Get order of hyperparameters
  if (is.null(order)) {
    hpar <- colnames(tuneGrid) # change order to change mapping
  } else {
    stopifnot(all(order %in% colnames(tuneGrid)))
    if (!all(colnames(tuneGrid) %in% order)) {
      order <- c(order, colnames(tuneGrid)[!colnames(tuneGrid) %in% order])
    }
    hpar <- order
  }
  
  ## Set plot mode
  switch(mode, 
         "size" = {
           cw <- case_when(ncol(tuneGrid) == 1 ~ {c("x",NA,NA,NA)},
                           ncol(tuneGrid) == 2 ~ {c("x","color",NA,NA)},
                           ncol(tuneGrid) == 3 ~ {c("x","color","size",NA)},
                           ncol(tuneGrid) >= 4 ~ {c("x","color","size","linetype")})
         },
         "alpha" = {
           cw <- case_when(ncol(tuneGrid) == 1 ~ {c("x",NA,NA,NA)},
                           ncol(tuneGrid) == 2 ~ {c("x","color",NA,NA)},
                           ncol(tuneGrid) == 3 ~ {c("x","color","alpha",NA)},
                           ncol(tuneGrid) >= 4 ~ {c("x","color","alpha","linetype")})
         })
  
  ## Plot performance metrics
  p.nrow <- ceiling(sqrt(n_distinct(met_train$.metric)))
  p <- 
    met_train %>%
    mutate(across(!!hpar[-c(which(cw == "x"), which(cw == "size!alpha"))], ~ factor(.x))) %>%
    ggplot(aes_string(x = hpar[1], y = "mean",
                      color = switch(!is.na(cw[2]),hpar[2],NULL),
                      linetype = switch(!is.na(cw[4]),hpar[4],NULL)))
  
  # add points
  if (is.na(cw[3])) {
    p <- p + 
      geom_point(size = point.size, alpha = point.alpha)
  } else {
    switch(mode,
           "size" = {
             p <- p +
               geom_point(aes_string(size = hpar[3]), alpha = point.alpha)
           },
           "alpha" = {
             p <- p +
               geom_point(aes_string(alpha = hpar[3]), size = point.size)
           })
  }
  
  # add lines
  if (is.na(cw[4])) {
    if (is.na(cw[3])) {
      gr <- switch(!is.na(cw[2]),hpar[2],NULL)
    } else {
      gr <- paste0("interaction(", paste0(c(hpar[2],hpar[3]), collapse =  ", "), ")")
    }
    p <- p +
      geom_line(aes_string(group = gr), alpha = line.alpha)
  } else {
    gr <- paste0("interaction(", paste0(c(hpar[2],hpar[3],hpar[4]), collapse =  ", "), ")")
    p <- p +
      geom_line(aes_string(group = gr), alpha = line.alpha)
  }
  if (y.log10) {
    p <- p + scale_y_log10(labels = scales::label_number())
  }
  if (x.log10) {
    if (is.numeric(met_train[[hpar[1]]])) {
      p <- p + scale_x_log10(labels = scales::label_number())
    } else {
      message(hpar[1], " not numeric, no log-transformation done!")
    }
  }
  p <- p + 
    geom_vline(aes(xintercept = bestTune[[hpar[1]]]), linetype = 2) +
    facet_wrap(~ .metric, scales = "fixed", nrow = p.nrow) +
    scale_color_viridis_d(option = "plasma", begin = .9, end = 0) +
    theme_bw() +
    ggtitle(paste0("Hyperparameter tuning: ", sub("_", "::", pkg_fun)))
  
  print(p)
  invisible(p)
}

print.tidytune <- function(tidytune) {
  
  stopifnot(inherits(tidytune, "tidytune"))
  
  cat("///////  tidytune  ///\n")
  mprec <- tidytune$mod$pre$actions$recipe$recipe
  mspec <- tidytune$mod$fit$actions$model$spec
  print(mprec)
  cat("\n\n")
  print(mspec)
  
  ntrain <- tidytune$nSplit %>% filter(set == "train") %>% dplyr::select("n") %>% as.numeric()
  ntest <- tidytune$nSplit %>% filter(set == "test") %>% dplyr::select("n") %>% as.numeric()
  ntot <- tidytune$nSplit %>% filter(set == "total") %>% dplyr::select("n") %>% as.numeric()
  
  ptrn <- paste0("# training samples: ", stringr::str_pad(ntrain, width = 9))
  ptst <- paste0("# test samples:     ", stringr::str_pad(ntest, width = 9))
  ptot <- paste0("# total samples:    ", stringr::str_pad(ntot, width = 9))
  ppar <- paste0("# parameters:       ", stringr::str_pad(tidytune$nX, width = 9))
  phyp <- paste0("# hyper-par config.:", stringr::str_pad(tidytune$nTune, width = 9))
  nhyp <- paste(names(tidytune$args$tuneGrid), collapse = " ; ")
  acc.test <-  tidytune$perf.test  %>% filter(.metric == "accuracy") %>% dplyr::select(".estimate")
  res <- list(spec = mspec,
              n = list(train = as.numeric(ntrain), test = ntest, total = ntot),
              p = list(nparam = as.numeric(tidytune$nX), nhyper = as.numeric(tidytune$nTune)),
              accuracy = as.numeric(acc.test))
  
  cat("> n                             > p\n")
  cat(paste(ptrn, ppar, sep = "   "), "\n")
  cat(paste(ptst, phyp, sep = "   "), "\n")
  cat(paste(ptot, nhyp, sep = "   "),"\n")
  
  cat("\n> performance\n")
  cat(paste0("# test set accuracy:", stringr::str_pad(acc.test, width = 9)), "\n")
  invisible(res)
}

summary.tidytune <- function(res) {
  
  stopifnot(inherits(res, "tidytune"))
  rec <- res$mod$pre$actions$recipe$recipe
  mod <- res$mod$fit$actions$model$spec
  
  y <- filter(rec$var_info, role == "outcome")$variable
  hpar <- colnames(tuneGrid)
  fform <- class(mod)[1]
  engine <- mod$engine
  transl <- translate(mod)
  pkg_fun <- paste(transl$method$fit$func, collapse = "_")
  
  metric.best <- res$args$metric.best
  mfun <- ifelse(attr(eval(parse(text = metric.best[1])),"direction") == "maximize", max, min)
  met_meta <- c(".metric",".estimator","mean","n","std_err",".config")
  
  tuned <- apply(res$args$tuneGrid, 2, function(x) {paste(range(x), collapse = ":")})
  tuned <- gsub(" ", "", paste0(paste(paste(names(tuned), tuned, sep = "["), collapse = "];"), "]"))
  
  best <- res$bestTune %>% dplyr::select(-any_of(met_meta))
  best <- paste0(paste(paste(names(best), best, sep = "["), collapse = "];"), "]")
  
  if (inherits(res$train, "tune_results")) {
    acc.train <- res$train %>% 
      show_best(metric = res$args$metric.best, n = 1) %>%
      dplyr::select(any_of(c(".metric","mean","n","std_err",".set")))
  } else {
    acc.train <- res$train$.metrics %>% 
      filter(.metric == res$args$metric.best) %>%
      filter(mean == mfun(mean)) %>%
      dplyr::select(any_of(c(".metric","mean","n","std_err",".set"))) %>%
      "["(1,)
  }
  
  acc.test <- res$perf.test %>% dplyr::select(c(".metric",".estimate"))
  acc.test.head <-t(acc.test)[1,] ; acc.test <- as.numeric(t(acc.test)[2,])
  acc.test <- array(data = acc.test, dim = c(1,length(acc.test)), dimnames = list(NULL, acc.test.head))
  acc.test <- as_tibble(acc.test) %>% rename_all(list(function(x) paste0("Test.", x)))
  
  dres <- 
    tibble(.rows = 1) %>%
    mutate(method = paste(res$fform, sub("_", "::", res$func)),
           tuned = tuned,
           best = best) %>%
    bind_cols(acc.train) %>%
    bind_cols(acc.test) %>%
    mutate(ntrain = res$nSplit %>% filter(set == "train") %>% dplyr::select("n") %>% as.numeric(),
           nvar = res$nX,
           nclass = rec$template %>% dplyr::select(all_of(y)) %>% unlist() %>% nlevels(),
           nhypercomb = unlist(res$nTune),
           elapsed = c(res$time$delta, res$time$elapsed),
           ncpu = res$time$ncpu)
  return(structure(dres, class = c("summary.tidytune", class(dres))))
}

predict.tidytune <- function(tidytune, new_data, ...) {
  
  stopifnot(inherits(tidytune, "tidytune"))
  
  # show_model_info(tidytune$fform) # returns NULL
  
  predict(tidytune$mod, new_data, ...)
  
}

# split data
split.data <- function(X, fac, group = NULL, by = NULL, SplitRatio = 2/3, nmin = 1,
                       drop = NULL, verbose = TRUE) {
  
  ## simon.crameri@usys.ethz.ch, Aug. 2021
  
  ## USAGE
  # X           data.frame      complete dataset to be split by rows
  # fac         factor          grouping factor of X
  # group       character       (optional) column name of X denoting a grouping variable. Rows of the same grouping variable cannot be split into different sets
  # by          character       (optional) column name of X denoting stratification variable. Rows of the same stratification variable tend to be split into different sets
  # SplitRatio  numeric         desired ratio of samples in training set
  # nmin        numeric         minimum number of rows per class (if group=NULL) or per group (if group!=NULL) in order to keep that class in training and test sets
  # drop        character       character string denoting unwanted classes (levels of <fac>), which will not appear in the training or test sets
  # verbose     logical         if TRUE, will print summary messages
  
  # check
  stopifnot(inherits(X, c("data.frame","matrix")),
            is.factor(fac), length(fac) == nrow(X),
            is.numeric(SplitRatio), SplitRatio >= 0, SplitRatio <= 1,
            is.numeric(nmin), nmin > 0)
  if (!is.null(group)) {
    stopifnot(is.character(group),
              group %in% colnames(X))
  }
  if (!is.null(by)) {
    stopifnot(is.character(by),
              by %in% colnames(X))
  }
  if (!is.null(drop)) {
    stopifnot(is.character(drop))
    if (!all(drop %in% levels(fac))) warning("These class(es) in <drop> are not in levels of <fac>:\n", paste(drop[!drop %in% levels(fac)], collapse = ","))
  }
  
  # remove small classes (n < nmin) or classes to be dropped
  if (!is.null(group)) {
    t <- sort(apply(table(fac, X[,group]), 1, function(x) {sum(x>0)}), decreasing = TRUE)
    names(t[t < nmin])
  } else {
    t <- table(fac)
    names(t[t < nmin])
  }
  fac2rm <- unique(c(names(t[t < nmin]), drop))
  X2 <- X[!fac %in% fac2rm,]
  fac2 <- droplevels(fac[!fac %in% fac2rm])
  dnew1 <- X[fac %in% fac2rm,]
  Ynew1 <- fac[fac %in% fac2rm]
  
  # split
  # sample = sample.split(Y = fac2, group = X2[,group], SplitRatio = SplitRatio) # caTools way
  sample = logical(length = nrow(X2))
  for (class in levels(fac2)) {
    
    # subset X by kept classes
    Xcl <- X2[fac2 %in% class,]
    
    # tabulate <group> and <by> factors
    if (!is.null(by)) {
      if (!is.null(group)) {
        t <- apply(table(Xcl[,group], Xcl[,by]), 2, as.numeric)
        rownames(t) <- names(table(Xcl[,group]))
      } else {
        t <- apply(table(rownames(Xcl), Xcl[,by]), 2, as.numeric)
        rownames(t) <- rownames(Xcl)
      }
    } else {
      if (!is.null(group)) {
        t <- apply(table(Xcl[,group], rep(1, nrow(Xcl))), 2, as.numeric)
        rownames(t) <- names(table(Xcl[,group]))
      } else {
        t <- apply(table(rownames(Xcl), rep(1, nrow(Xcl))), 2, as.numeric)
        rownames(t) <- rownames(Xcl)
      }
    }
    
    # select most frequent <group> per <by> factor, add to training set proportionally to SplitRatio
    i <- apply(t, 2, function(x) {
      xord <- x[order(x, decreasing = T)]
      xmax <- xord[xord == max(xord)]
      x1 <- sample(names(xmax), 1)
      
      xrest <- x[!names(x) %in% x1]
      ct <- xrest[xrest>0]
      
      r <- (SplitRatio*length(x[x>0])) - 1
      x2 <- names(ct[sample(seq(length(ct)), round(r))])
      
      c(x1, x2)
    })
    
    # remove some training samples if the effective proportion is too high
    sel <- as.vector(unlist(i))
    
    if (length(sel) > round(nrow(t) * SplitRatio)) {
      # cat(class)
      d <- names(which.min(rowSums(t[i[[which.max(lengths(i))]],,drop=F])))
      sel <- sel[!sel %in% d]
    }
    
    # update sample vector
    if (!is.null(group)) {
      w <- Xcl[,group] %in% sel
    } else {
      w <- rownames(Xcl) %in% sel
    }
    sample[fac2 %in% class] <- w
  } 
  
  # create training and validation sets
  dtrain = X2[sample,]
  dtest  = X2[!sample,]
  
  Ytrain <- droplevels(fac2[sample])
  Ytest <- droplevels(fac2[!sample])
  
  # # correct training and validation sets (only needed if caTools::sample.split() is used)
  # train = X2[sample,]
  # test  = X2[!sample,]
  # switch(method,
  #        "test2train" = {
  #          dtest <- test[!test[,group] %in% train[,group],]
  #          dtrain <- X2[!rownames(X2) %in% rownames(dtest),]
  # 
  #          Ytest <- droplevels(fac2[!sample][!test[,group] %in% train[,group]])
  #          Ytrain <- droplevels(fac2[!rownames(X2) %in% rownames(dtest)])
  #        },
  #        "train2test" = {
  #          dtrain <- train[!train[,group] %in% test[,group],]
  #          dtest <- X2[!rownames(X2) %in% rownames(dtrain),]
  # 
  #          Ytrain <- droplevels(fac2[sample][!train[,group] %in% test[,group]])
  #          Ytest <- droplevels(fac2[!rownames(X2) %in% rownames(dtrain)])
  #        })
  
  # create new set (not in training or validation set)
  dnew <- dtest[!Ytest %in% Ytrain,]
  Ynew <- factor(Ytest[!Ytest %in% Ytrain], levels = levels(fac))
  
  dtest <- dtest[Ytest %in% Ytrain,]
  Ytest <- droplevels(Ytest[Ytest %in% Ytrain])
  
  dnew <- rbind(dnew1, dnew)
  Ynew <- factor(c(as.character(Ynew1), as.character(Ynew)))
  
  # check
  if (!is.null(group)) {
    stopifnot(!any(dtrain[,group] %in% dtest[,group]))
  }
  stopifnot(nrow(dtrain) == length(Ytrain),
            nrow(dtest) == length(Ytest),
            all(Ytrain %in% Ytest),
            all.equal(levels(Ytrain), levels(Ytest)))
  
  # make summary
  if ("fac" %in% colnames(X)) {
    dtrain[,"fac.replaced"] <- dtrain[,"fac"]
    dtest[,"fac.replaced"] <- dtest[,"fac"]
    dnew[,"fac.replaced"] <- dnew[,"fac"]
  }
  dtrain$fac <- Ytrain
  dtest$fac <- Ytest
  dnew$fac <- Ynew
  dp <- data.frame(ntrain = as.numeric(table(Ytrain)),
                   ntest = as.numeric(table(Ytest)),
                   grouptrain = as.numeric(sapply(levels(Ytrain), function(x) {length(unique(dtrain[dtrain[,"fac"] %in% x,group]))})),
                   grouptest = as.numeric(sapply(levels(Ytest), function(x) {length(unique(dtest[dtest[,"fac"] %in% x,group]))})),
                   ptrain = as.matrix(table(Ytrain) / table(droplevels(fac2[fac2 %in% Ytrain]))),
                   ptest = as.matrix(table(Ytest) / table(droplevels(fac2[fac2 %in% Ytrain]))))
  
  # be verbose
  if (verbose) {
    message("Training set: ", nrow(dtrain), "; Test set: ", nrow(dtest), "; New set: ", nrow(dnew))
    message("\nDesired SplitRatio: ", round(SplitRatio, 2), " ; Effective SplitRatio: ", round(nrow(dtrain)/(nrow(dtrain)+nrow(dtest)), 2), " (", paste(round(range(dp$ptrain), 2), collapse = " - "), ")")
    message("\nDropped these ", length(table(Ynew)), " / ", length(table(fac)),
            " classes from training and test set:\n", 
            paste(sort(names(table(Ynew))), collapse = ","))
    message("\nKept these ", length(table(Ytrain)), " / ", length(table(fac)),
            " classes in training and test set:\n",
            paste(sort(names(table(Ytrain))), collapse = ","))
  }
  
  # return
  args <- list(group = group, by = by, SplitRatio = SplitRatio, nmin = nmin, drop = drop)
  invisible(list(Xtrain = dtrain, Xtest = dtest, Xnew = dnew,
                 Ytrain = Ytrain, Ytest = Ytest, Ynew = Ynew,
                 args = args, summary = dp))
}
funky <- function(n) {
  ramp <- grDevices::colorRampPalette(c("#A6CEE3","#1F78B4","#B2DF8A",
                                        "#33A02C","#FB9A99","#E31A1C",
                                        "#FDBF6F","#FF7F00","#CAB2D6",
                                        "#6A3D9A","#FFFF99","#B15928"))
  ramp(n)
}
funky2 <- function(n) {
  ramp <- grDevices::colorRampPalette(c("#A6CEE3","#1F78B4","#B2DF8A",
                                        "#33A02C","#FB9A99","#E31A1C",
                                        "#FDBF6F","#FF7F00","#CAB2D6",
                                        "#6A3D9A","#FFFF99","#B15928"))
  c(ramp(n-1), "gray30")
}
spectre <- function(n) {
  ramp <- grDevices::colorRampPalette(c("#A6CEE3","#438EC0","#63A8A0","#98D277",
                                        "#E6E600","#B89B74","#F16667","#7A0403",
                                        "#B5B5B5","#FE982C","#3BA432","#C3AAD2",
                                        "#7D54A5","#00204D","#EAD27A", "#B15928"))
  ramp(n)
}
spectre2 <- function(n) {
  ramp <- grDevices::colorRampPalette(c("#A6CEE3","#438EC0","#63A8A0","#98D277",
                                        "#E6E600","#B89B74","#F16667","#7A0403",
                                        "#B5B5B5","#FE982C","#3BA432","#C3AAD2",
                                        "#7D54A5","#00204D","#EAD27A", "#B15928"))
  c(ramp(n-1), "gray30")
}
# n <- 16
# pie(rep(1,n),col=spectre(n))
# pie(rep(1,n+1),col=spectre2(n+1))

get.leaflet <- function(sp, map.layers = NULL, idfac = "ID_Lab", classfac = "Species", layerfac = NULL, popup = c(idfac, classfac, layerfac),
                        link = NULL, radius = 4, opacity = 0.6, jitter.amount = 0.001, print = TRUE,
                        providers = c("OpenStreetMap.Mapnik","Esri.WorldImagery","Esri.WorldTopoMap","Esri.WorldGrayCanvas"),
                        palette = adegenet::funky(nlevels(sp@data[,classfac]))) {
  library(leaflet)
  
  ## Arguments
  # sp  spatialPointsDatFrame
  # map.layers  list of data map layers. each element is a character vector of the form c(<polygon>,<variable>,<palette>). The polygon and palette will be accessed using get()
  # idfac character denoting point label variable 
  # classfac character denoting point color variable
  # layerfac character denoting overlayGroups factor variable
  # popup character vector denoting point variables used for popups
  # link character denoting variable with SpecimenID
  # palette character vector with point colours
  
  ## Check input
  stopifnot(inherits(sp, "SpatialPointsDataFrame"),
            idfac %in% colnames(sp@data),
            classfac %in% colnames(sp@data),
            all(popup %in% colnames(sp@data)))
  if (!is.factor(sp@data[,classfac])) {sp@data[,classfac] <- factor(sp@data[,classfac])}
  if (!is.null(layerfac)) {
    stopifnot(layerfac %in% colnames(sp@data))
    if (!is.factor(sp@data[,layerfac])) sp@data[,layerfac] <- factor(sp@data[,layerfac])
  }
  if (!is.null(link)) {
    stopifnot(link %in% colnames(sp@data))
    linkbase = "https://tropicos.org/Specimen/"
    sp@data["LINK"] <- paste0("<a href='", linkbase, sp@data[,link], "', target=\"blank\">", sp@data[,link], "</a>")
    popup <- unique(c("LINK", popup))
  }
  
  ## Helperfunctions
  get.popup <- function(sp, popup) {
    unname(gsub("^</br>", "", apply(
      matrix(sapply(seq(length(popup)), FUN = function(x) {
        apply(sp@data[,popup[x],drop=F], MARGIN = 1, FUN = function(y) {paste0(paste0("</br>", popup[x], ": "), y)})
      }), ncol = length(popup)), MARGIN = 1, FUN = paste, collapse = "")))
  }
  
  ## Add tiles
  M <- leaflet()
  for (i in 1:length(providers)) {
    M <- M %>% 
      addProviderTiles(providers[i], group = providers[i])
  }
  
  ## Jitter points
  set.seed(4065)
  sp@coords <- apply(sp@coords, 2, jitter, factor = 0, amount = jitter.amount)
  
  ## Add layers
  map.groups <- character()
  for (i in map.layers) {
    M <- M %>% 
      addPolygons(data = get(i[1]),
                  group = i[1],
                  popup = get(i[1])[[i[2]]],
                  fillColor = get(i[3])(get(i[1])[,i[2]][[1]]))
    map.groups <- c(map.groups, i[1])
  }
  
  ## Add jittered oints
  pal <- colorFactor(palette = palette, domain = sp@data[,classfac])
  if (is.null(layerfac)) {
    M <- M %>% 
      addCircleMarkers(data = sp, color = pal(sp@data[,classfac]),
                       radius = radius, opacity = opacity, fillOpacity = opacity, stroke = FALSE,
                       popup = get.popup(sp, popup),
                       label = sp@data[,idfac],
                       group = classfac)
  } else {
    for (i in 1:nlevels(droplevels(sp@data[,layerfac]))) {
      layer <- levels(sp@data[,layerfac])[i]
      sp.layer <- sp[sp@data[,layerfac] %in% layer,]
      M <- M %>%
        addCircleMarkers(data = sp.layer, color = pal(sp.layer@data[,classfac]),
                         radius = radius, opacity = opacity, fillOpacity = opacity, stroke = FALSE,
                         popup = get.popup(sp.layer, popup),
                         label = sp.layer@data[,idfac],
                         group = layer)
    }
  }
  
  ## Add layer control
  if (is.null(layerfac)) {
    M <- M %>% addLayersControl(baseGroups = providers,
                                overlayGroups = c(map.groups, classfac),
                                position = "topleft")
  } else {
    M <- M %>% addLayersControl(baseGroups = providers,
                                overlayGroups = c(map.groups, levels(sp@data[,layerfac])),
                                position = "topleft") %>%
      hideGroup(levels(sp@data[,layerfac])[-1])
  }
  if (print) print(M)
  invisible(M)
}

## Transform
eco.transform <- function(df, oldnames, newnames,
                          var.geol = "ECO_surfaceLithologyAfrica",
                          var.vege = "ECO_vegetationAtlas2007",
                          var.dist = c("ECO_dist2Coast","ECO_dist2InlandWaters"),
                          var.pa = "ECO_ProtectedAreas", var.forest = "ECO_forestCover2017",
                          var.elev = "ECO_elevation", elev_div100 = TRUE,
                          temp_div10 = TRUE) {
  
  # rename
  names(df)[match(oldnames, names(df))] <- newnames
  
  # # geol
  # legend.geol <- c("0 unmapped",
  #                  "1 Carbonate [Marble (Cipolin)]",
  #                  "2 Karst [Mesozoic Limestones incl Tsingy]",
  #                  "3 NonCarbonate [Sandstones]",
  #                  "4 Metasedimentary [Basement Rocks]",
  #                  "5 Alkaline Intrusive Volcanic",
  #                  "6 Silicic [Basement Rocks]",
  #                  "7 Metaigneous [Basement Rocks]",
  #                  "8 Ultramafic [Ultrabasics]",
  #                  "9 Extrusive Volcanic [Lavas incl Basalts & Gabbros]",
  #                  "11 Hydric Organic [Mangroves]",
  #                  "18 Alluvium [Fluvial, Beach, Unconsolidated Sands]", # includes 14, 15
  #                  "20 Water")
  # 
  # dgeol <- data.frame(geol = sapply(strsplit(legend.geol, split = " "), "[", 1),
  #                     legend.geol = legend.geol)
  # 
  # df[,var.geol] <- factor(dgeol[match(as.character(df[,var.geol]), dgeol$geol),"legend.geol"])
   
  # geology: combine classes
  g <- df[,var.geol]
  comb1 <- levels(g)[grep(" Alluvium ", levels(g))] # combine 14,15,18
  repl1 <- "18 Alluvium [Various]"
  comb2 <- levels(g)[grep("Metasedimentary|Silicic|Metaigneous", levels(g))] # combine 4,6,7
  repl2 <- "4 Basement Rocks [Metasedimentary, Silicic, Metaigneous]"
  new <- replace(replace(as.character(g), g %in% comb1, repl1), g %in% comb2, repl2)
  newlev <- levels(factor(new))[order(as.numeric(sapply(strsplit(levels(factor(new)), split = " "), "[", 1)))]
  df[,var.geol] <- factor(new, levels = newlev)
  
  # # vege
  # legend.vege <- c("0 unmapped / Clouds",
  #                  "1 Water Bodies",
  #                  "2 Bare soil/rock",
  #                  "3 Mangroves",
  #                  "4 Cultivation",
  #                  "5 Western dry forest [NW/W]",
  #                  "6 Plateau grassland-wooded grassland mosaic",
  #                  "7 Wooded grassland-bushland mosaic",
  #                  "9 Western humid forest",
  #                  "10 Western dry forest [SW]",
  #                  "11 Degraded south western dry spiny forest",
  #                  "12 South western dry spiny forest-thicket",
  #                  "13 Wetlands",
  #                  "14 Humid forest",
  #                  "15 Littoral forest",
  #                  "16 Degraded humid Forest",
  #                  "18 unmapped / Clouds",
  #                  "19 South western coastal bushland",
  #                  "22 Western sub-humid forest",
  #                  "23 Tapia forest",
  #                  "25 Sea")
  # 
  # dvege <- data.frame(vege = sapply(strsplit(legend.vege, split = " "), "[", 1),
  #                     legend.vege = legend.vege)
  # 
  # df[,var.vege] <- factor(dvege[match(as.character(df[,var.vege]), dvege$vege),"legend.vege"])
  
  # vege: combine classes
  v <- df[,var.vege]
  comb1 <- levels(v)[grep("Clouds", levels(v))]
  repl1 <- "0 Unmapped / Clouds"
  comb2 <- levels(v)[grep("Western dry forest", levels(v))]
  repl2 <- "5 Western dry forest"
  new <- replace(replace(as.character(v), v %in% comb1, repl1), v %in% comb2, repl2)
  newlev <- levels(factor(new))[order(as.numeric(sapply(strsplit(levels(factor(new)), split = " "), "[", 1)))]
  df[,var.vege] <- factor(new, levels = newlev)
  
  # Bio01 to Bio11: divide by 10 degrees if temp_div10
  for (i in newnames[grep("temp|isothermality", newnames)]) {
    if (temp_div10) {
      df[,i] <- round(df[,i] / 10, 2)
    } else {
      df[,i] <- round(df[,i], 2)
    }
  }

  # Bio12 to Bio19 unchanged
  for (i in newnames[grep("prec", newnames)]) {
    df[,i] <- round(df[,i], 2)
  }
  
  # elevation divide by 100
  if (elev_div100) {
    df[,var.elev] <- round(df[,var.elev]/100)
  } else {
    df[,var.elev] <- round(df[,var.elev])
  }
  
  # Distance unchanged
  for (i in var.dist) {
    df[,i] <- round(df[,i], 0)
  }
  
  # NA to zero: forest
  df[,var.forest] <- as.numeric(as.character(df[,var.forest]))
  df[,var.forest][is.na(df[,var.forest])] <- 0
  df[,var.forest] <- as.factor(df[,var.forest])
  
  # NA to zero: Protected Areas  
  df[,var.pa] <- as.character(df[,var.pa])
  df[,var.pa][is.na(df[,var.pa])] <- "unprotected"
  df[,var.pa] <- factor(df[,var.pa], levels = c("I-IV","V-VI","treat_as_V-VI","MARINE","unprotected"))
  
  # return results
  return(df)
}
recode.wa <- function(df, split = "", rm.redundant = FALSE, na2zero = FALSE, cor.thr = 0.999) {
  ncol1 <- ncol(df)
  fac <- sapply(df, is.character) | sapply(df, is.factor)
  make_dummy <- names(fac[fac])
  df <- cbind(dummy(df[,make_dummy], split = split, rm.redundant = rm.redundant, na2zero = na2zero),
              df[,!names(df) %in% make_dummy])
  ncol2 <- ncol(df)
  message("Recoded ", ncol1, " variables to ", ncol2, " variables (including dummy)")
  
  # remove redundant variables (one out of two binary dummy variables)
  dcor <- suppressWarnings(cor(df, method = "spearman", use = "everything"))
  diag(dcor) <- dcor[upper.tri(dcor)] <- NA
  lcor <- apply(dcor, 2, function(x) {which(abs(x) >= cor.thr)})
  redundant <- names(which(lengths(lcor) > 0))
  message(paste0("Found ", length(redundant), " / ", ncol(df), " redundant variables: ", paste(redundant, collapse = ",")))
  
  # return results
  res <- list(df = df, redundant = redundant)
  return(res)
}
annotate.meta <- function(dx, dy, fac, qual, min.frac.dummy = 0.05) {
  
  stopifnot(all.equal(rownames(dx), rownames(dy)))
  
  ## Usage:
  # dx              data.frame qualitative and quantitative data
  # dy              data.frame metadata
  # fac             FCT grouping factor
  # qual            CHR vector of names of qualitative characters that need to be converted to dummy variables
  # min.frac.dummy  NUM only create dummy variables for character states that occur with at least this frequency
  
  # bind
  dy <- cbind(dx, dy[,!names(dy) %in% names(dx)])
  
  # annotate
  dy[,"Species_N"] <- factor(paste0("D. ", sub("([A-Za-z. ])_n([0-9]+)", "\\1 (n = \\2", annotate.n(fac)), ")"))
  dy[,"ZT_Dalbergia_Supergroup"] <- ifelse(dy[,"ZT_Dalbergia_Group"] %in% c("Chapelieri","Maritima","Pervillei","Tricolor"), "SGI",
                                    ifelse(dy[,"ZT_Dalbergia_Group"] %in% "Bracteolata", "Bracteolata",
                                    ifelse(dy[,"ZT_Dalbergia_Group"] %in% "Xerophila", "Xerophila", "SGII")))
  
  # factorize
  for (i in c("ZT_Dalbergia_Supergroup","ZT_Dalbergia_Group","Species_N")) {dy[,i] <- factor(dy[,i])}
  
  # make dummy
  for (i in qual) {
    if (i %in% names(dy)) {
      tab.dummy <- sort(table(dy[,i], useNA = "always"), decreasing = TRUE)
      id.dummy <- names(tab.dummy[tab.dummy > sum(tab.dummy) * min.frac.dummy])
      dy <- factor2dummy(df = dy, factorname = i, keep.levels = id.dummy)
    }
  }
  return(dy)
}
df2pca <- function(df, fac = NULL) {
  
  stopifnot(all(sapply(df, is.numeric)))
  
  # imput missing values
  if (!is.null(fac)) {
    dd.pca <- data.frame(apply(df, MARGIN = 2, FUN = impute.na, method = "group.mean", groups = fac))
  } else {
    dd.pca <- df
  }
  dd.pca <- data.frame(apply(dd.pca, MARGIN = 2, FUN = impute.na, method = "mean"))
  
  ## Remove monomorphic variables
  mono <- names(which(apply(dd.pca, 2, function(x) var(x) == 0)))
  n1 <- ncol(dd.pca)
  message("found ", length(mono), " / ", n1, " monomorphic variables: ", paste(mono, collapse = ","))
  dd.pca <- dd.pca[,!names(dd.pca) %in% c(mono)]
  message("keeping ", ncol(dd.pca), " / ", n1, " variables for multivariate analyses")
  return(dd.pca)
}
recode.IAWA <- function(strings, names, split = "-", code_variable = "v$", code_unsure = "\\?$",
                        lookup, lookup.idx = NULL, value = 1, value_variable = 0.5, value_unsure = NA,
                        var.qual = NULL, var.quant = NULL) {
  
  ## Usage
  # strings   CHR  vector of IAWA strings
  # names     CHR  vector of names (same length as <strings>)
  # split     CHR  character used to split <strings>
  # code_variable CHR regex used to find variable character states
  # code_unsure   CHR regex used to find unsure character states
  # lookup    DF   lookup data.frame used to convert IAWA strings to recoded characters
  # lookup.idx NUM vector of 3 numbers, indicating column number in <lookup> matching IAWA, RECODED, and QUANTITATIVE_VALUE_ASSIGNED
  # value     NUM  value assigned if character state is present and not variable
  # value_variable NUM value assigned if character state is present and variable
  # var.qual  CHR  character vector denoting characters (must be in lookup[,lookup.idx[2]])
  # var.quant CHR  character vector denoting quantitative characters (must be in lookup[,lookup.idx[2]])
  
  ## Author: simon.crameri@usys.ethz.ch, Dec 2021
  
  # guess 
  if (is.null(lookup.idx)) {
    guess1 <- sapply(lookup, is.numeric)
    guess2 <- apply(lookup, 2, function(x) any(grepl("_",x)) & all(nchar(na.omit(x)) < 8))
    lookup.idx <- c(which(guess1)[1], which(guess2)[1])
    guess3 <- sapply(lookup[!grepl("_", lookup[,lookup.idx[2]]) & !is.na(lookup[,lookup.idx[2]]),], function(x) {
      y <- suppressWarnings(as.numeric(x)) ; length(na.omit(y)) > 0
    })
    lookup.idx <- c(lookup.idx, which(guess3)[!which(guess3) %in% lookup.idx][1])
  }
  if (is.null(var.qual)) {
    vars <- as.character(na.omit(lookup[,lookup.idx[2]]))
    var.qual <- unique(vars[grepl("_", vars)])
  } else {
    var.qual <- unique(var.qual)
  }
  if (is.null(var.quant)) {
    vars <- as.character(na.omit(lookup[,lookup.idx[2]]))
    var.quant <- unique(vars[!grepl("_", vars)])
  } else {
    var.quant <- unique(var.quant)
  }
  
  # check input
  stopifnot(grepl(split, strings),
            length(names) == length(strings),
            inherits(lookup, c("matrix","data.frame")),
            is.numeric(value),
            is.numeric(value_variable) | is.na(value_variable),
            is.numeric(value_unsure) | is.na(value_unsure),
            is.numeric(lookup[,lookup.idx[1]]),
            length(suppressWarnings(na.omit(as.numeric(lookup[,lookup.idx[3]])))) > 0,
            any(var.qual %in% lookup[,lookup.idx[2]]),
            any(var.quant %in% lookup[,lookup.idx[2]]))
                       
  
  # helperfunction
  .recode.IAWA <- function(string, split = "-", code_variable = "v$", code_unsure = "?$", lookup) {
    ss <- unlist(strsplit(string, split = "-"))
    sv <- grepl(code_variable, ss)
    su <- grepl(code_unsure, ss)
    ss <- gsub(code_unsure, "", gsub(code_variable, "", ss))
    sr <- sapply(as.numeric(ss), FUN = function(x) {lookup[lookup[,lookup.idx[1]] == x,lookup.idx[2]]})
    return(data.frame(code = ss, code_variable = sv, code_unsure = su, recoded = sr))
  }
  
  # list character table for each string
  ll <- lapply(as.character(strings), FUN = function(x) {
    .recode.IAWA(string = x, split = split, code_variable = code_variable,
                 code_unsure = code_unsure, lookup = lookup)
  })

  # convert to data.frame  
  df <- data.frame(array(data = NA, dim = c(length(ll), length(c(var.qual,var.quant))),
                            dimnames = list(NULL, c(var.qual,var.quant))))
  
  # fill qualitative and quantitative
  for (i in seq(length(ll))) {
    t <- ll[[i]]
    for (j in var.qual) {
      if (j %in% t[,"recoded"]) {
        if (all(t[t[,"recoded"] %in% j,"code_variable"])) {
          df[i,j] <- value_variable
        } else {
          if (all(t[t[,"recoded"] %in% j,"code_unsure"])) {
            df[i,j] <- value_unsure
          } else {
            df[i,j] <- value
          }
        }
      }
    }
    for (j in var.quant) {
      if (j %in% t[,"recoded"]) {
        dcodes <- t[t[,"recoded"] %in% j,]
        fills <- numeric()
        for (code in dcodes$code) {
          # r <- gsub("^-", "", gsub("[A-Za-z µ]", "", lookup[lookup[,1] %in% code,3]))
          # fill <- suppressWarnings(c(fill, mean(na.omit(as.numeric(unlist(strsplit(r, split = "-|>|<|=")))))))
          if (dcodes[dcodes[,"code"] == code,"code_unsure"]) {
            fill <- value_unsure
          } else {
            fill <- as.numeric(lookup[lookup[,lookup.idx[1]] %in% code,lookup.idx[3]])
          }
          fills <- c(fills, fill)
        }
        df[i,j] <- mean(fills)
      }
    }
  }
  
  # fill zeros (assume that a character was assessed if one or more character states were coded, fill absent states of assessed character with zero)
  var.qual2 <- unique(gsub("_[a-z]", "", var.qual))
  for (i in var.qual2) {
    vars <- var.qual[grep(paste0("^", i, "_"), var.qual)]
    rec <- t(apply(df[,vars,drop=F], MARGIN = 1, FUN = function(x) {if (anyNA(x) & sum(x,na.rm=T)>0) x[is.na(x)] <- 0 ; x}))
    if (nrow(rec) == nrow(df)) df[,vars] <- rec else df[,vars] <- t(rec)
  }
  
  # name
  if (any(duplicated(names))) {
    repeat{
      names[duplicated(names)] <- paste0(names[duplicated(names)], "_R")
      if (!any(duplicated(names))) break()
    }
  }
  names(ll) <- names
  rownames(df) <- names
  
  # return everything
  args <- list(strings = strings, names = names, split = split, code_variable = code_variable,
               lookup = lookup, lookup.idx = lookup.idx, value = value, value_variable = value_variable,
               value_unsure = value_unsure, var.qual = var.qual, var.quant = var.quant)
  res <- list(list = ll, df = df, args = args)
  return(res)
}
factor2dummy <- function(df, factorname, keep.levels = levels(factor(df[,factorname])), na.method = "zero") {
  
  dfac <- data.frame(factor(df[,factorname]))
  names(dfac) <- factorname
  keep.levels <- levels(dfac[,factorname])[levels(dfac[,factorname]) %in% keep.levels]
  
  m <- model.matrix(~dfac[,factorname])
  m <- m[match(rownames(dfac), rownames(m)), , drop = FALSE]
  m[,1] <- ifelse(dfac[,factorname] == levels(dfac[,factorname])[1], 1, 0) # first column is a dummy intercept (first level not returned in model.matrix)
  colnames(m) <- paste0(factorname, "_", levels(dfac[,factorname]))
  rownames(m) <- rownames(dfac)
  if (na.method == "zero") m[is.na(m)] <- 0
  
  stopifnot(all.equal(sum(m), sum(table(dfac[,factorname]))))
  df2 <- df[,!colnames(df) %in% factorname]
  pos <- which(colnames(df) == factorname)
  if (pos == 1) {
    df2 <- cbind(m[,paste0(factorname, "_", keep.levels),drop=FALSE], df2)
  } else if (pos == ncol(df)) {
    df2 <- cbind(df2, m[,paste0(factorname, "_", keep.levels),drop=FALSE])
  } else {
    df2 <- cbind(df2[,1:(pos-1),drop=FALSE], m[,paste0(factorname, "_", keep.levels),drop=FALSE], df2[,pos:ncol(df2),drop=FALSE])
  }
  names(df2)[pos:(pos + length(keep.levels) - 1)] <- paste0(factorname, "_", keep.levels)
  rownames(df2) <- rownames(df)
  return(df2)
}
dummy <- function(.data, split = FALSE, count = FALSE, rm.redundant = TRUE, na2zero = FALSE) {
  if (!"dplyr" %in% installed.packages()) install.packages("dplyr")
  if (is.null(dim(.data))) .data <- data.frame(dummy = .data)
  .dummy <- function(x, split = FALSE, count = FALSE, rm.redundant = TRUE, na2zero = FALSE) {
    
    # identify factor levels and combinations
    l <- strsplit(unlist(as.character(x)), split = split)
    lev <- sort(unique(unlist(l)))
    if (count) {
      z <- lapply(l, FUN = function(y) {table(factor(y, levels = lev))})
    } else {
      z <- lapply(l, FUN = function(y) {lev %in% y})
    }
    
    # handle NA's
    if (count) {
      if (!na2zero & anyNA(x)) {
        z[[which(is.na(x))]] <- setNames(rep(NA, length(lev)), lev)
      }
    } else {
      if (!na2zero & any(is.na(l))) {
        # z[[which(is.na(l))]] <- rep(NA, length(lev))
        z[is.na(l)] <- lapply(z[is.na(l)], function(x) rep(NA, length(x)))
      }
    }
    
    # create dummy variables
    if (count) {
      m <- dplyr::bind_rows(z)
    } else {
      m <- matrix(as.numeric(unlist(z)), nrow = length(unlist(x)), byrow = TRUE)
    }
    colnames(m) <- paste0("_", lev)
    
    # remove redundant dummy variables
    if (rm.redundant) {
      if (count) {
        redundant <- character()
      } else {
        dcor <- cor(m)
        diag(dcor) <- dcor[upper.tri(dcor)] <- NA
        lcor <- apply(dcor, 2, function(x) {which(x == -1)})
        redundant <- names(which(lengths(lcor) > 0))
      }
      m <- m[,!colnames(m) %in% redundant, drop = FALSE]
    }
    
    # return results
    return(m)
  }
  M <- .data[,0]
  for (i in seq(ncol(.data))) {
    m <- .dummy(x = .data[,i], split = split, count = count,
                rm.redundant = rm.redundant, na2zero = na2zero)
    colnames(m) <- paste0(colnames(.data)[i], colnames(m))
    M <- dplyr::bind_cols(M, as.data.frame(m), .name_repair = "minimal")
  }
  return(M)
}

split.fun <- function(x, labs, digits, varlen, faclen) {
  labs <- gsub("<", "\n<", gsub(">", "\n>", gsub(" ", "\n", gsub(":", "\n", labs))))
  for (i in 1:length(labs)) {
    
    # split labs[i] into multiple lines
    labs[i] <- paste(strwrap(labs[i], width = 10), collapse = "\n")
    rm(i)
  }
  labs
}
get.diag <- function(df, fac, level = nlevels(fac)-1, verbose = TRUE) {
  
  ## Usage:
  # df      data.frame
  # fac     grouping factor (same length as nrow(df))
  # level   grouping level. 1 only searches diagnostic features in a single group, nlevels(fac)-1 [default] considers all possible combinations of grouping factor levels.
  # verbose whether to be verbose during comupations [default: TRUE]
  
  ## Value:
  # data.frame with the following components
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
  
  ## Author: simon.crameri@usys.ethz.ch, Jan 2022
  
  # check input
  if (!is.factor(fac)) fac <- factor(fac)
  if (level >= nlevels(fac)) level <- nlevels(fac)-1
  stopifnot(inherits(df, c("matrix","data.frame")),
            all.equal(length(fac), nrow(df)),
            is.numeric(level), level >=1)
  
  # search qualitative and quantitative variables
  bin <- suppressWarnings(sapply(names(df), function(x) is.logical(df[,x]) | (is.numeric(df[,x]) & all(na.omit(df[,x]) <=1) & all(na.omit(df[,x]) >= 0)) | (is.factor(df[,x]) & length(unique(x)) == 2)))
  v.qual <- names(bin[bin])
  v.quant <- names(bin[!bin])
  
  # be verbose
  fname <- gsub("__", "_", gsub(" ", "_", gsub('^factor|[()]|[]]|[[]|"|`', "", sub("(.*\\\\)(.*)", "(\\2)", sub("(.*\\$)(.*)", "(\\2)", gsub(",", "_", deparse(substitute(fac))))))))
  if (verbose) {
    message("Searching diagnostic characters to distinguish ", fname, "\n  ", length(v.qual), " qualitative and ", length(v.quant), " quantitative variables at level ", level, ".")
  }
  
  # loop over variables using different frequency thresholds
  res <- data.frame(FACTOR = fname, NLEVELS = nlevels(fac), LEVEL = NA, GROUP = NA,
                    VARIABLE = NA, THR = NA, STATE = NA, FREQ = NA, MISS = NA, TPR = NA, FPR = NA)[0,]
  
  # loop over variables, value thresholds (quantitative only) and frequency thresholds
  for (var in c(v.qual, v.quant))  {
    
    v <- df[,var]
    
    # check
    if (is.factor(v)) v <- as.numeric(v)
    stopifnot(!all(is.na(v)), is.numeric(v))
    if (var %in% v.qual) stopifnot(na.omit(v) <= 1, na.omit(v >= 0))
    
    # set thresholds for quantitative variables
    if (var %in% v.quant) {
      r <- range(v[!v %in% min(v, na.rm = T)], na.rm = T)
      thrs.range <- seq(r[1], r[2], length.out = 1000) # if more than 1000 unique values, use range
      thrs.uniq <- sort(unique(v))[-1]
      thrs <- list(thrs.range, thrs.uniq)[[which.min(c(length(thrs.range), thrs.uniq))]]
    } else {
      thrs <- 1 
    }
    
    # loop over value thresholds (quantitative only) and frequency thresholds 
    for (thr in thrs) {
      
      # tabulate
      if (var %in% v.qual) {
        t <- table(Group = fac, Variable = v)
      } else {
        t <- table(Group = fac, Variable = factor(v >= thr, levels = c(FALSE, TRUE)))
        colnames(t) <- c("0","1")
      }
      
      # skip if variable is constant
      if (!all(c("0","1") %in% colnames(t))) {
        # if (verbose) cat("skipping invariant", var, "\r")
        next()
      } else {
        # if (verbose) cat(var,"\r")
      }
      
      for (freq in seq(1, 0.5 , by = -0.01)) {
        
        # high frequency or always present in
        high <- apply(t, 1, function(x) {(x[length(x)]/sum(x)) >= freq})
        
        # low frequency or always absent in
        low <- apply(t, 1, function(x) {(x[1]/sum(x)) >= freq})
        
        l <- length(low)
        shigh <- sum(high)
        slow <- sum(low)
        
        fhigh <- shigh/l
        flow <- slow/l
        
        # test for diagnosticity at level 1 (single group), 2 (two groups), 3 (three groups), etc.
        for (lev in seq(level)) {
          lowlev <- flow > 0 & flow <= lev/l
          highlev <- fhigh > 0 & fhigh <= lev/l
          assign(x = paste0("low", lev), value = lowlev)
          assign(x = paste0("high", lev), value = highlev)
          rm(list = c("lev","lowlev","highlev"))
        }
        
        # calculate TPR and FPR
        f <- fac[!is.na(v)]
        m <- sum(is.na(v))/nrow(df)
        
        clow <- names(low[which(low)])
        tpr.low <- sum(t[levels(fac) %in% clow,"0"])/sum(f %in% clow)
        fpr.low <- sum(t[!levels(fac) %in% clow,"0"])/sum(!f %in% clow)
        
        chigh <- names(high[which(high)])
        tpr.high <- sum(t[levels(fac) %in% chigh,"1"])/sum(f %in% chigh)
        fpr.high <- sum(t[!levels(fac) %in% chigh,"1"])/sum(!f %in% chigh)
        
        # save result if diagnostic at level <level>
        if (get(paste0("low", level))) {
          res.low <- data.frame(FACTOR = fname, NLEVELS = nlevels(fac), LEVEL = length(clow),
                                GROUP = paste(clow, collapse = "+"),
                                VARIABLE = var, THR = thr, STATE = "<", FREQ = freq, MISS = m,
                                TPR = tpr.low, FPR = fpr.low)
          res <- rbind(res, res.low)
        }
        if (get(paste0("high", level))) {
          res.high <- data.frame(FACTOR = fname, NLEVELS = nlevels(fac), LEVEL = length(chigh), 
                                 GROUP = paste(chigh, collapse = "+"),
                                 VARIABLE = var, THR = thr, STATE = "≥", FREQ = freq, MISS = m,
                                 TPR = tpr.high, FPR = fpr.high)
          res <- rbind(res, res.high)
        }
      }
    }
  }
  
  # order results
  res <- res[order(res$FACTOR, res$VARIABLE, res$GROUP, -res$TPR, res$FPR),]
  
  # filter reduntant results
  res <- res[!duplicated(as.matrix(res[,!names(res) %in% "FREQ"])),]
  res <- res[!duplicated(as.matrix(res[,!names(res) %in% c("THR","FREQ","TPR","FPR")])),]
  
  # be verbose
  if (verbose) {
    ndiag <- length(unique(subset(res, TPR >= 0.9 & FPR <= 0.1)$VARIABLE))
    ndiaggrp <- length(unique(subset(res, TPR == 1 & FPR == 0 & LEVEL == 1)$VARIABLE))
    message("  found ", ndiag, " potentially useful (TPR >= 90%, FPR <= 10%) variables,\n  of which ", ndiaggrp, " are diagnostic (TPR = 100%, FPR = 0%) to distinguish single groups!")
  }
  
  # return results
  rownames(res) <- NULL
  res$THR[res$VARIABLE %in% v.qual] <- NA
  res$STATE[res$VARIABLE %in% v.qual] <- ifelse(res$STATE[res$VARIABLE %in% v.qual] == "<", "absent", "present")
  return(res)
}
# correlation
varCorr <- function(cor, cor.thr = 0.9, dev = NULL) {
  
  # check
  stopifnot(length(unique(dim(cor))) == 1,
            !anyNA(cor[upper.tri(cor)]))
  # if (!is.null(names(dev))) {
  #   stopifnot(all.equal(names(dev)[1:ncol(cor)], colnames(cor)))
  # }
  
  # keep upper triangle only
  cor[lower.tri(cor)] <- diag(cor) <- NA
  
  # define groups of correlated varialbes
  l.cor <- apply(cor, 1, function(x) {w <- which(abs(x) > cor.thr) ; if (length(w) > 0) colnames(cor)[w] else NA})
  
  # loop through groups and find variables to be eliminated (keep variable with high deviance gain)
  correlated <- character()
  for (i in seq(length(l.cor))) {
    
    # variables
    corr.base <- names(l.cor)[i]
    corr.with <- l.cor[[i]]
    
    if (!is.null(dev)) {
      # criterion: deviances (good criterion, invariant to scale of predictor)
      dev.base <- dev[pmatch(corr.base, names(dev))] # dev[i]
      dev.with <- dev[pmatch(corr.with, names(dev))]
      dev.max <- c(corr.base, corr.with)[which.max(c(dev.base, dev.with))]
      torm <- c(corr.base, corr.with)[!c(corr.base, corr.with) %in% dev.max]
    } else {
      # criterion: use the first variable in the correlation matrix order
      torm <- corr.with[!corr.with %in% names(l.cor)[1:i]]
    }    
    correlated <- as.character(na.omit(unique(c(correlated, torm))))
  }
  return(correlated)
}
varFilter <- function(X, cor.thr = 0.9, cor.method = "spearman", rank = NULL,
                      freqCut = 95/5, uniqueCut = 10, verbose = T) {
  
  # find variables with near-zero variation (monomorphic)
  nzv <- caret::nearZeroVar(x = X, freqCut = freqCut, uniqueCut = uniqueCut)
  monomorphic <- colnames(X)[nzv]
  
  # find variables with high collinearity
  dcor <- cor(X[,colnames(X)[!colnames(X) %in% monomorphic]], method = cor.method)
  if (is.null(rank)) dev <- rank else dev <- -1 * rank
  if (!is.null(dev) & any(!colnames(X) %in% names(dev))) {
    miss <- colnames(dcor)[!colnames(X) %in% names(dev)]
    warning("found ", length(miss), " variables missing in <rank>, these will be discarded: ", paste(miss, collapse = ","))
  }
  correlated <- varCorr(cor = dcor, cor.thr = cor.thr, dev = dev)
  
  # filter out monomorphic / highlX correlated variables
  torm <- union(monomorphic, correlated)
  message(paste0("found ", length(monomorphic), "/", ncol(X), " (", round(100*length(monomorphic)/ncol(X),2), "%) highly monomorphic variables"))
  if (verbose) print(monomorphic)
  message(paste0("found ", length(correlated), "/", ncol(X), " (", round(100*length(correlated)/ncol(X),2), "%) highly correlated variables"))
  if (verbose) print(correlated)
  message(paste0("filtered ", length(torm), "/", ncol(X), " (", round(100*length(torm)/ncol(X),2), "%) variables"))
  v <- names(X)[!names(X) %in% torm]
  
  # select variables for downstream analysis
  message(paste0("selected ", length(v), "/", ncol(X), " (", round(100*length(v)/ncol(X), 2), "%) variables"))
  if (verbose) print(v)
  invisible(X[,v])
}
subset.cor <- function(df, method = "spearman", use = "pairwise", thr = 0.7) {
  dc <- cor(df, method = method, use = use)
  diag(dc) <- NA
  dc2 <- dc
  dc[upper.tri(dc)] <- NA
  dc2[lower.tri(dc2)] <- NA
  s <- apply(dc, MARGIN = 2, function(x) {any(abs(x[!is.na(x)]) >= thr)})
  s2 <- apply(dc2, MARGIN = 2, function(x) {any(abs(x[!is.na(x)]) >= thr)})
  dc[s|s2,s|s2]
}