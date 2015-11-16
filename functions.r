### Larry Fenn
### larry.fenn@gmail.com
### 08/14/2015
### Parallel R code for multivariate Bayesian linear regression on Windows in R
###
### Check commands.r for more information about these functions.

reqPkgs <- c("snow", "parallel", "BLR")
new.packages <- reqPkgs[!(reqPkgs %in% installed.packages()[,"Package"])]
if (length(new.packages)) {install.packages(new.packages)}

library(snow)
library(parallel)
library(BLR)

coreCount <- detectCores() ## Parallel machine count. Change this at your peril!
centralizeCap <- 1000 ## Number of models to centralize.
centralizeRepeat <- 7 ## Number of repeat iterations for centralization.

### pre:  x is a data frame where each column is an independent variable, each
###       over the same number of rows.
###       level is a correlation significance level.
###       dimensions is the maximum number of independent variables used.
### post: A matrix where each column is viable "model"; stored as indices of the
###       independent variables.
modelSelect <- function(x, level = .6, dimensions) {
    if (dimensions == 1) {
        return(matrix(1:ncol(x), 1))
    }
    varCombos <- combn(1:ncol(x), dimensions)

    ## Parallel machine initialization.
    cl <- makeCluster(coreCount, type = "SOCK")
    clusterExport(cl, "level", env = environment())

    ## Evaluate on clusters a function which returns TRUE if the correlation
    ## matrix between all variables in the combination has entries all (in
    ## absolute value) below "level", except the diagonal (since a variable is
    ## 100% correlated with itself).
    results <- parLapplyLB(cl, 1:ncol(varCombos), function(i) {
            length(which(abs(cor(x[,varCombos[,i]])) >= level)) == dimensions
        })

    ## Close parallel machines.
    stopCluster(cl)

    ## The return type of parLapplyLB is a list object where each entry is a
    ## result; hence, unlist() must be used to flatten the list into a vector.
    varCombos[,which(unlist(results))]
}

### pre:  Two data frames, hopefully nonempty, that we wish to glue together by
###       columns (i.e. the number of columns in the result is the sum of the
###       number of columns in the input), attaching an appropriate amount of
###       NA entries to the ends to fill in the gaps.
### post: The joined data frame. the ordering of columns is preserved.
joinFrames <- function(x, y) {
    if (is.null(nrow(x)) || ncol(x) == 0) {
        return(y)
    }
    if (is.null(nrow(y)) || ncol(y) == 0) {
        return(x)
    }
    reversed <- FALSE
    if (nrow(x) < nrow(y)) {
        z <- x
        x <- y
        y <- z
        reversed <- TRUE
    }
    y <- rbind(y, matrix(NA, nrow(x) - nrow(y), ncol(y)))
    if (reversed) {
        return(cbind(y, x))
    }
    cbind(x, y)
}

### pre:  A list of models (each row is a model, whose entries are indices
###       corresponding to the independent variable in indvar in the model),
###       the independent variable list itself, and the dependent variable
###       time series.
### post: An output list of the following format:
###       Entry 1: BLR results data frame.
###       Entry 2: BLR results data frame, sorted by average squared residuals.
###       Entry 3: Centralized BLR results of the top models.
###       Entry 4: Centralized BLR results, sorted by average squared residuals.
###       Entry 5: Run time on BLR step (excluding centralization).
###       Entry 6: Run time on the whole process.
generateResults <- function(modelList, indvar, depvar) {
    startTime <- proc.time()
    print("Now computing BLR...")

    ## Initialize parallel: send each cluster the namespace and packages.
    cl <- makeCluster(4, type = "SOCK")
    ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
    clusterExport(cl, ex)
    clusterExport(
        cl,
        c("modelList", "indvar", "depvar"),
        env = environment()
        )

    clusterEvalQ(cl, library(BLR))

    ## Evaluate on clusters the helper function and output into a data frame.
    results <- parLapplyLB(cl, 1:nrow(modelList), blrHelper)
    ## The return type of parLapplyLB is a list object where each entry is a
    ## result; hence, unlist() must be used to flatten the list into a vector.
    results <- data.frame(matrix(unlist(results), ncol = 35, byrow = T))

    results[,17] <- rank(results[,14], ties.method = "min")
    results[,18] <- rank(results[,15], ties.method = "min")
    ## Column 19 is ranked in descending order.
    results[,19] <- rank(
        -as.numeric(levels(results[,16])[results[,16]]),
        ties.method = "min"
        )
                             ## Column Name             Index
    colnames(results) <- cbind("Variable_Lag",       ## 1
                               "Variable_1",         ## 2
                               "Variable_2",         ## 3
                               "Variable_3",         ## 4
                               "Variable_4",         ## 5
                               "Variable_5",         ## 6
                               "Intercept",          ## 7
                               "Coefficient_Lag",    ## 8
                               "Coefficient_1",      ## 9
                               "Coefficient_2",      ## 10
                               "Coefficient_3",      ## 11
                               "Coefficient_4",      ## 12
                               "Coefficient_5",      ## 13
                               "AvgSqResiduals",     ## 14
                               "DevianceIC",         ## 15
                               "PseudoRSquared",     ## 16
                               "Ranked_ASR",         ## 17
                               "Ranked_DIC",         ## 18
                               "Ranked_PRS",         ## 19
                               "lm_Intercept",       ## 20
                               "lm_Coefficient_Lag", ## 21
                               "lm_Coefficient_1",   ## 22
                               "lm_Coefficient_2",   ## 23
                               "lm_Coefficient_3",   ## 24
                               "lm_Coefficient_4",   ## 25
                               "lm_Coefficient_5",   ## 26
                               "lm_Pvalue_Intercept",## 27
                               "lm_Pvalue_Lag",      ## 28
                               "lm_Pvalue_1",        ## 29
                               "lm_Pvalue_2",        ## 30
                               "lm_Pvalue_3",        ## 31
                               "lm_Pvalue_4",        ## 32
                               "lm_Pvalue_5",        ## 33
                               "lm_Rsquared",        ## 34
                               "lm_AdjRsquared"      ## 35
                               )

    timingPhase1 <- proc.time() - startTime
    print(timingPhase1)
    print("Now centralizing...")

    sortOrder <- order(results[,17]) #Avg. sq. residuals
    resultsSorted <- results[sortOrder,] ## Avg. sq. residuals.
    numCentralize <- min(nrow(resultsSorted), centralizeCap)
    centInput <- modelList[sortOrder,]

    ## Evaluate on clusters the helper function and output into a data frame.
    clusterExport(cl, "centralizeRepeat")
    centralizeResults <- parLapplyLB(
        cl,
        1:numCentralize,
        function(i) {centralizeHelper(centInput[i,])}
        )
    centralizeResults <- data.frame(
        matrix(unlist(centralizeResults), ncol = 35, byrow = T)
        )

    stopCluster(cl) ## Retire the clusters.

    centralizeResults[,17] <- rank(centralizeResults[,14], ties.method = "min")
    centralizeResults[,18] <- rank(centralizeResults[,15], ties.method = "min")
    centralizeResults[,19] <- rank(
        -as.numeric(levels(centralizeResults[,16])[centralizeResults[,16]]),
        ties.method = "min"
        )

    colnames(centralizeResults) <- colnames(results)

    ## The remaining columns are inherited (unmodified by centralizing).
    for (j in 20:35) {
        centralizeResults[,j] <- resultsSorted[1:numCentralize,j]
    }
    centralizeResultsSorted <- centralizeResults[order(centralizeResults[,17]),]

    timingPhase2 <- proc.time() - startTime
    print(timingPhase2)
    print("Finished.")

    output <- list(
        "OptResults" = results,
        "OptResultsSorted" = resultsSorted,
        "OptCentralizeResults" = centralizeResults,
        "OptCentralizeResultsSorted" = centralizeResultsSorted,
        "Runtime" = timingPhase1,
        "Runtime2" = timingPhase2
        )
    return(output)
}

### pre:  The environment variable 'modelList' exists, with a list of models.
###       Additionally, indvar and depvar exist as environment variables (these
###       should have been sent to the environments with clusterExport).
### post: A vector output encoding the BLR output.
blrHelper <- function(modelIndex) {
    ## Grab all of the current model independent variables, up to 5 (if needed).
    modelFactors <- cbind(indvar[,modelList[modelIndex,1]])
    if (!is.na(modelList[modelIndex,2])) {
        modelFactors <- cbind(modelFactors, indvar[,modelList[modelIndex,2]])
    }
    if (!is.na(modelList[modelIndex,3])) {
        modelFactors <- cbind(modelFactors, indvar[,modelList[modelIndex,3]])
    }
    if (!is.na(modelList[modelIndex,4])) {
        modelFactors <- cbind(modelFactors, indvar[,modelList[modelIndex,4]])
    }
    if (!is.na(modelList[modelIndex,5])) {
        modelFactors <- cbind(modelFactors, indvar[,modelList[modelIndex,5]])
    }

    sink("NUL") ## Suppress the BLR output to console.
    varB <- BLR(depvar, modelFactors)
    sink()
    sqResiduals <- (varB$yHat - varB$y)^2
    avgSqResidual <- mean(subset(sqResiduals, !is.na(sqResiduals)))

    pseudoR <- cor(varB$yHat, varB$y)^2  ## Efron's pseudo-R-squared.
    numVars <- 5 - sum(is.na(modelList[modelIndex,]))
    ## Adjust pseudo-R-squared for number of variables.
    pseudoRadj  <- 1 - ((1 - pseudoR) * ((length(depvar) - 1) /
                                             (length(depvar) - numVars - 1)))

    lmTemp <- lm(depvar ~ modelFactors)
    lmTemp <- summary(lmTemp)

    ## Result packaging step.
    result <- rep(NA, 35)

    result[1]  <- NA ## First lag: not used at the moment.
    ## Model variable names.
    if (!is.na(modelList[modelIndex,1])) {
        result[2]  <- names(indvar)[modelList[modelIndex,1]]
    }
    if (!is.na(modelList[modelIndex,2])) {
        result[3]  <- names(indvar)[modelList[modelIndex,2]]
    }
    if (!is.na(modelList[modelIndex,3])) {
        result[4]  <- names(indvar)[modelList[modelIndex,3]]
    }
    if (!is.na(modelList[modelIndex,4])) {
        result[5]  <- names(indvar)[modelList[modelIndex,4]]
    }
    if (!is.na(modelList[modelIndex,5])) {
        result[6]  <- names(indvar)[modelList[modelIndex,5]]
    }
    result[7]  <- varB$mu ## Intercept.
    result[8]  <- NA ## Coefficient on the first lag (not used).
    ## Model variable coefficients.
    if (!is.na(modelList[modelIndex,1])) {
        result[9]  <- varB$bF[1]
    }
    if (!is.na(modelList[modelIndex,2])) {
        result[10] <- varB$bF[2]
    }
    if (!is.na(modelList[modelIndex,3])) {
        result[11] <- varB$bF[3]
    }
    if (!is.na(modelList[modelIndex,4])) {
        result[12] <- varB$bF[4]
    }
    if (!is.na(modelList[modelIndex,5])) {
        result[13] <- varB$bF[5]
    }
    result[14] <- avgSqResidual
    result[15] <- varB$fit$DIC
    result[16] <- pseudoRadj
    result[20] <- coef(lmTemp)[1,1] ## OLS intercept.
    result[21] <- NA ## OLS first lag coefficient (not used).
    ## OLS variable coefficients.
    if (!is.na(modelList[modelIndex,1])) {
        result[22] <- coef(lmTemp)[2,1]
    }
    if (!is.na(modelList[modelIndex,2])) {
        result[23] <- coef(lmTemp)[3,1]
    }
    if (!is.na(modelList[modelIndex,3])) {
        result[24] <- coef(lmTemp)[4,1]
    }
    if (!is.na(modelList[modelIndex,4])) {
        result[25] <- coef(lmTemp)[5,1]
    }
    if (!is.na(modelList[modelIndex,5])) {
        result[26] <- coef(lmTemp)[6,1]
    }
    result[27] <- coef(lmTemp)[1,4] ## OLS P-value on the intercept.
    result[28] <- NA ## OLS P-value on the first lag (not used).
    ## OLS variable P-values.
    if (!is.na(modelList[modelIndex,1])) {
        result[29] <- coef(lmTemp)[2,4]
    }
    if (!is.na(modelList[modelIndex,2])) {
        result[30] <- coef(lmTemp)[3,4]
    }
    if (!is.na(modelList[modelIndex,3])) {
        result[31] <- coef(lmTemp)[4,4]
    }
    if (!is.na(modelList[modelIndex,4])) {
        result[32] <- coef(lmTemp)[5,4]
    }
    if (!is.na(modelList[modelIndex,5])) {
        result[33] <- coef(lmTemp)[6,4]
    }
    result[34] <- lmTemp$r.squared
    result[35] <- lmTemp$adj.r.squared

    result
}

### pre:  The environment variable 'modelList' exists, with a list of models.
###       varIndices is given as the indices of the variables to be centralized.
###       Additionally, indvar and depvar exist as environment variables (these
###       should have been sent to the environments with clusterExport).
### post: A vector output encoding the BLR output of centralized models.
centralizeHelper <- function(varIndices) {
    ## Grab all of the current model independent variables, up to 5 (if needed).
    modelFactors <- cbind(indvar[,varIndices[1]])
    if (!is.na(varIndices[2])) {
        modelFactors <- cbind(modelFactors, indvar[,varIndices[2]])
    }
    if (!is.na(varIndices[3])) {
        modelFactors <- cbind(modelFactors, indvar[,varIndices[3]])
    }
    if (!is.na(varIndices[4])) {
        modelFactors <- cbind(modelFactors, indvar[,varIndices[4]])
    }
    if (!is.na(varIndices[5])) {
        modelFactors <- cbind(modelFactors, indvar[,varIndices[5]])
    }

    coefLog <- matrix(NA, centralizeRepeat, 6)
    dicLog <- matrix(NA, centralizeRepeat, 1)

    for (j in 1:centralizeRepeat) {
        sink("NUL") ## Suppress BLR output to console.
        varB <- BLR(depvar, modelFactors)
        sink()
        coefLog[j,1] <- varB$mu
        coefLog[j,2] <- varB$bF[1]
        coefLog[j,3] <- varB$bF[2]
        coefLog[j,4] <- varB$bF[3]
        coefLog[j,5] <- varB$bF[4]
        coefLog[j,6] <- varB$bF[5]
        dicLog[j,1]  <- varB$fit$DIC
    }

    modelCoeffs <- c(
        median(coefLog[,2]),
        median(coefLog[,3]),
        median(coefLog[,4]),
        median(coefLog[,5]),
        median(coefLog[,6])
        )

    ## Only use as many coefficients as there are model factors (strip NAs).
    modelCoeffs <- as.vector(matrix(
        modelCoeffs,
        length(modelCoeffs) - sum(is.na(modelCoeffs)),
        1
        ))
    yhatNoIntercept <- t(modelCoeffs * t(modelFactors))
    yhat <- as.matrix(rowSums(yhatNoIntercept)) + median(coefLog[,1])
    sqResiduals <- (yhat - depvar) ^ 2
    avgSqResiduals <- mean(sqResiduals)

    pseudoR <- cor(yhat, depvar) ^ 2 ## Efron's pseudo R-squared.
    numVar <- 5 - sum(is.na(varIndices))
    ## Adjust for number of variables.
    pseudoRadj  <- 1 - ((1 - pseudoR) * ((length(depvar) - 1) /
                                             (length(depvar) - numVar - 1)))

    ## Result packaging step.
    result <- rep(NA, 35)

    result[1]  <- NA ## First lag: not used at the moment.
    ## Model variable names inherited from the input.
    result[2]  <- names(indvar)[varIndices[1]]
    result[3]  <- names(indvar)[varIndices[2]]
    result[4]  <- names(indvar)[varIndices[3]]
    result[5]  <- names(indvar)[varIndices[4]]
    result[6]  <- names(indvar)[varIndices[5]]
    ## Model coefficients and intercept from centralization.
    result[7]  <- median(coefLog[,1])
    result[8]  <- NA ## First lag not used.
    result[9]  <- median(coefLog[,2])
    result[10] <- median(coefLog[,3])
    result[11] <- median(coefLog[,4])
    result[12] <- median(coefLog[,5])
    result[13] <- median(coefLog[,6])
    result[14] <- avgSqResiduals
    result[15] <- median(dicLog[,1])
    result[16] <- pseudoRadj

    result
}

### pre:  A list for output of the format that generateResults creates, and
###       a name that will be used for the output directory.
### post: Output .csv files generated in the ./name/ directory.
outputBLRresults <- function(x, name) {
    dir.create(paste("./", name, sep = ""), showWarnings = FALSE)
    for (i in 1:4) {
        write.csv(
            x[i],
            paste("./", name, "/", names(x)[i], ".csv", sep = ""),
            row.names = FALSE
            )
    }
    write.csv(x$Runtime[3], paste("./", name, "/", "Runtime.csv", sep = ""))
    write.csv(x$Runtime2[3], paste("./", name, "/", "Runtime2.csv", sep = ""))
}
