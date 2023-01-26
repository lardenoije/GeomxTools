####################################################################################################
# Run a linear model on a GeoMxSet                                                                 #
# author: Roy Lardenoije                                                                           #
# last updated: 20-07-2022                                                                         #
####################################################################################################

# Function based on the mixedModelDE function from the GeomxTools package. That function only allows
# for mixed models, which is not necessary when doing comparisons within a single slide. The current
# function works basically the same, but uses lm instead of lmer and thus accepts models without
# random effects.

modelDE <- function (object, elt = "exprs", modelFormula = NULL, groupVar = "group", 
          nCores = 1, multiCore = TRUE, pAdjust = "BY") 
{
    if (is.null(modelFormula)) {
        modelFormula <- design(object)
    }
    mTerms <- all.vars(modelFormula)
    if ("1" %in% mTerms) {
        mTerms <- mTerms[which(!(mTerms %in% "1"))]
    }
    if (!groupVar %in% mTerms) {
        stop("Error: groupVar needs to be defined as fixed effect in the model.\n")
    }
    if (any(!mTerms %in% names(sData(object)))) {
        stop("Error: Not all terms in the model formula are in pheno or protocol data.\n")
    }
    pDat <- sData(object)[, mTerms, drop=FALSE]
    for (i in names(pDat)) {
        if (inherits(i, "character")) {
            pDat[, i] <- as.factor(pDat[, i])
        }
    }
    if (nCores > 1) {
        deFunc <- function(i, groupVar, pDat, modelFormula, exprs) {
            dat <- data.frame(expr = exprs$exprs[i, ], pDat)
            lmOut <- suppressWarnings(summary(lm(modelFormula, dat))$coefficients)
            lmOut <- matrix(cbind(lmOut[-1, "Estimate", drop=FALSE], lmOut[-1, "Pr(>|t|)", drop=FALSE]), 
                             ncol = 2, dimnames = list(gsub(groupVar, "", rownames(lmOut)[-1]), 
                                                       c("Estimate", "Pr(>|t|)")))
            return(lmOut)
        }
        exprs <- new.env()
        exprs$exprs <- assayDataElement(object, elt = elt)
        if (multiCore & Sys.info()["sysname"] != "Windows") {
            lmOut <- parallel::mclapply(featureNames(object),
                                        deFunc, groupVar, pDat, 
                                        formula(paste("expr", as.character(modelFormula)[2], sep = " ~ ")),
                                        exprs, mc.cores = nCores)
        }
        else {
            cl <- parallel::makeCluster(getOption("cl.cores", nCores))
            lmOut <- parallel::parLapply(cl, featureNames(object), deFunc, groupVar, pDat, 
                                         formula(paste("expr", as.character(modelFormula)[2], sep = " ~ ")),
                                         exprs)
            suppressWarnings(parallel::stopCluster(cl))
        }
        lmOut <- do.call("rbind", lmOut)
        rownames(lmOut) <- featureNames(object)
    }
    else {
        deFunc <- function(expr, groupVar, pDat, modelFormula) {
            dat <- data.frame(expr = expr, pDat)
            lmOut <- suppressWarnings(summary(lm(modelFormula, dat))$coefficients)
            lmOut <- matrix(cbind(lmOut[-1, "Estimate", drop=FALSE], lmOut[-1, "Pr(>|t|)", drop=FALSE]), 
                            ncol = 2, dimnames = list(gsub(groupVar, "", rownames(lmOut)[-1]), 
                                                      c("Estimate", "Pr(>|t|)")))
            return(lmOut)
        }
        lmOut <- t(assayDataApply(object, 1, deFunc, groupVar, pDat, 
                                formula(paste("expr", as.character(modelFormula)[2], sep = " ~ ")), 
                                elt = elt))
    }
    if (!is.null(pAdjust)) {
        lmOut <- cbind(lmOut, p.adjust(lmOut[, 2], method = pAdjust))
    }
    colnames(lmOut) <- c("Estimate", "Pr(>|t|)", pAdjust)
    return(lmOut)
}