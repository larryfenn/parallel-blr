### Larry Fenn
### larry.fenn@gmail.com
### 08/14/2015
### Parallel R code for multivariate Bayesian linear regression on Windows in R.
###
### Directions: Run the script. The user will be prompted for further input.
###
### commands.r: The script that takes the user input and outputs into
###             appropriate directories the results of the analysis.
###             The user will be asked for two file inputs, and results will be
###             stored in folders (one for each dependent variable).
###             Data format: each column must have a header, and there should
###             not be a time column (i.e. the first column is the first data
###             series). There must be the same number of rows in all files.
###
### functions.r: The functions used by commands.r.
###     + coreCount, centralizeCap, and centralizeRepeat are hardcoded values
###       that can be changed in functions.r to suit whatever needs exist.
###     - modelSelect: Returns all possible lists of five independent variables
###       where pairwise correlation between the variables is below a given
###       threshold level (variable modelCorLevel).
###     - joinFrames: Helper function to join two data frames over fixed rows
###     - generateResults: First, runs Bayesian linear regression over all
###       possible models. Next, the top models are selected and the median of
###       multiple Bayesian regressions is taken as the final result.
###     - blrHelper: Helper function that is evaluated on parallel machines.
###     - centralizeHelper: Helper that is evaluated on parallel machines.
###     - ouptutBLRresults: I/O function to output  generateResults.
###
### TODO: Currently does not have first lag support.
### Note: First-time runs through the code may pop up Windows Firewall messages.
###       It does not matter if R is allowed access to the network for the
###       script to work.

source("functions.r")

modelCorLevel <- .1 ## Correlation level for model selection.

## Alternatively, hardcode the filenames here.
print("Select independent variables source:")
indvar <- read.csv(file.choose())
print("Select dependent variables source:")
depvarlist <- read.csv(file.choose())

modelList <- numeric(0)
print("Selecting models...")

for (i in 1:5) {
    modelList <- joinFrames(modelList, modelSelect(indvar, modelCorLevel, i))
}
modelList <- t(modelList)

if (ncol(modelList) < 5) {
    modelList <- cbind(modelList, matrix(
        NA,
        nrow(modelList),
        5 - ncol(modelList)
        ))
}

print(paste(nrow(modelList), "models under consideration."))

for (i in 1:ncol(depvarlist)) {
    print(paste("Now working on ", colnames(depvarlist)[i], "...", sep = ""))
    output <- generateResults(modelList, indvar, unlist(depvarlist[i]))
    outputBLRresults(output, colnames(depvarlist)[i])
}
print("Done!")
