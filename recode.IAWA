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
