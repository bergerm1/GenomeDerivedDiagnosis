library(testthat)
library(MolecularDiagnosis)
library(data.table)
context("GDD tests")

test_that("generate_feature_table", {
  example_data_dir <- system.file(package = "MolecularDiagnosis", "example_data")
  zipfile <- file.path(example_data_dir, "msk_impact_2017.zip")
  data_dir <- file.path(example_data_dir, "msk_impact_2017")
  feature_table_filename <- file.path(example_data_dir, "feature_table.tsv")

  utils::unzip(zipfile = zipfile, exdir = example_data_dir)

  arg_line <- paste(data_dir, feature_table_filename)
  generate_feature_table(arg_line)

  unlink(data_dir, recursive=TRUE)

  expect_true(file.exists(feature_table_filename))
})

test_that("train_classifier", {
  example_data_dir <- system.file(package = "MolecularDiagnosis", "example_data")
  feature_table_filename <- file.path(example_data_dir, "feature_table.tsv")
  model_filename <- file.path(example_data_dir, "test_classifier.rds")

  arg_line <- paste(feature_table_filename, model_filename, "--test")
  train_classifier(arg_line)
  expect_true(file.exists(model_filename))
})

test_that("predict_new", {
  example_data_dir <- system.file(package = "MolecularDiagnosis", "example_data")
  model_filename <- file.path(example_data_dir, "test_classifier.rds")
  feature_table_filename <- file.path(example_data_dir, "feature_table.tsv")
  prediction_table_filename <- file.path(example_data_dir, "prediction_table.xlsx")

  arg_line <- paste(model_filename, feature_table_filename, prediction_table_filename)
  predict_new(arg_line)
  expect_true(file.exists(prediction_table_filename))

  prediction_table <- as.data.table(readxl::read_excel(prediction_table_filename))

  expect_equal(dim(prediction_table), c(10945L, 51L))
  expect_equal(
    head(prediction_table$SAMPLE_ID),
    c("P-0000004-T01-IM3", "P-0000015-T01-IM3", "P-0000023-T01-IM3",
      "P-0000024-T01-IM3", "P-0000025-T01-IM3", "P-0000025-T02-IM5")
  )
  expect_equal(
    head(prediction_table$Conf1),
    c(0.151, 0.4409, 0.3189, 0.168, 0.1651, 0.1932),
    tolerance=1e-4
  )

  expect_equal(prediction_table[Classification_Category == "train", mean(Cancer_Type == Pred1)], 0.3075343, tolerance=1e-4)
})

## start with default rfcal

test_that("rfcal default options", {

  expect_match(rfcal$label, "Random Forest with Probability Calibration")
  expect_match(rfcal$library[1], "randomForest")
  expect_match(rfcal$library[2], "glmnet")
  expect_equal(rfcal$loop, NULL)
  expect_equal(length(rfcal$type), 2)
  expect_equal(length(rfcal$parameters$parameter), 4)
  expect_equal(length(rfcal$parameters$class), 4)
  expect_equal(length(rfcal$parameters$label), 4)
  ## grid()
  ## fit()
  ## predict()
  ## prob()
  ## predictors()
  ## varImp()
  ## levels()
  ## tags()
  ## sort()
  ## oob()

})


## tdt
test_that("check tdt() ", {
  ## input
  iris_dt = as.data.table(iris)
  expect_equal(dim(iris_dt)[1], 150)
  expect_equal(dim(iris_dt)[2], 5)
  expect_equal(length(colnames(iris_dt)), 5)
  ## check dimensions
  expect_equal(dim(tdt(iris_dt))[1], 4)
  expect_equal(dim(tdt(iris_dt))[2], 151)
  expect_equal(length(tdt(iris_dt)[,1]$Sepal.Length), 4)
})



### rename_types


### get_fc


### predict_new



## test GGDmetrics, based on multiClassSummary
## https://github.com/topepo/caret/blob/master/pkg/caret/tests/testthat/test_multiclassSummary.R


test_that("test GDDmetrics ", {

  N = 1000
  M = 2
  set.seed(42)
  xTrain = matrix(runif(N*M), nrow = N)
  colnames(xTrain) = sapply(1:M, function(u){paste0(collapse = '', letters[sample(26, 3, replace = TRUE)])})

  yTrain = as.factor(letters[sample(3, N, replace = TRUE)])

  trCntlListMulti = trainControl(method = "cv", number = 3, verboseIter = FALSE, classProbs = TRUE, summaryFunction = GDDmetrics)

  enFitMulti = train(x = xTrain, y = yTrain, trControl = trCntlListMulti, method = "knn", tuneLength = 2)

  ## trCntlListMulti
  expect_equal(length(names(trCntlListMulti)), 27)
  ##  [1] "method"            "number"            "repeats"
  ##  [4] "search"            "p"                 "initialWindow"
  ##  [7] "horizon"           "fixedWindow"       "skip"
  ## [10] "verboseIter"       "returnData"        "returnResamp"
  ## [13] "savePredictions"   "classProbs"        "summaryFunction"
  ## [16] "selectionFunction" "preProcOptions"    "sampling"
  ## [19] "index"             "indexOut"          "indexFinal"
  ## [22] "timingSamps"       "predictionBounds"  "seeds"
  ## [25] "adaptive"          "trim"              "allowParallel"

  expect_match(trCntlListMulti$method, "cv")
  expect_equal(trCntlListMulti$number, 3)
  expect_equal(trCntlListMulti$repeats, NA)
  expect_match(trCntlListMulti$search, "grid")
  expect_equal(trCntlListMulti$p, 0.75)
  expect_equal(trCntlListMulti$initialWindow, NULL)
  expect_equal(trCntlListMulti$horizon, 1)
  expect_true(trCntlListMulti$fixedWindow)
  expect_equal(trCntlListMulti$skip, 0)
  expect_false(trCntlListMulti$verboseIter)
  expect_true(trCntlListMulti$returnData)
  expect_match(trCntlListMulti$returnResamp, "final")
  expect_false(trCntlListMulti$savePredictions)
  expect_true(trCntlListMulti$classProbs)
  ## trCntlListMulti$summaryFunction
  expect_match(trCntlListMulti$selectionFunction, "best")
  ## trCntlListMulti$preProcOptions
  ## [1] "thresh"    "ICAcomp"   "k"         "freqCut"   "uniqueCut" "cutoff"
  expect_equal(trCntlListMulti$preProcOptions$thresh, 0.95)
  expect_equal(trCntlListMulti$preProcOptions$ICAcomp, 3)
  expect_equal(trCntlListMulti$preProcOptions$k, 5)
  expect_equal(trCntlListMulti$preProcOptions$freqCut, 19)
  expect_equal(trCntlListMulti$preProcOptions$uniqueCut, 10)
  expect_equal(trCntlListMulti$preProcOptions$cutoff, 0.9)
  expect_equal(trCntlListMulti$sampling, NULL)
  expect_equal(trCntlListMulti$index, NULL)
  expect_equal(trCntlListMulti$indexOut, NULL)
  expect_equal(trCntlListMulti$indexFinal, NULL)
  expect_equal(trCntlListMulti$timingSamps, 0)
  expect_equal(trCntlListMulti$predictionBounds, c(FALSE, FALSE))
  expect_equal(trCntlListMulti$seeds, NA)
  ##trCntlListMulti$adaptive
  ## [1] "min"      "alpha"    "method"   "complete"
  expect_equal(trCntlListMulti$adaptive$min, 5)
  expect_equal(trCntlListMulti$adaptive$alpha, 0.05)
  expect_match(trCntlListMulti$adaptive$method, "gls")
  expect_true(trCntlListMulti$adaptive$complete)
  expect_false(trCntlListMulti$trim)
  expect_true(trCntlListMulti$allowParallel)

  ## enFitMulti
  expect_equal(length(names(enFitMulti)), 20)
  ## [1] "method"       "modelInfo"    "modelType"    "results"      "pred"
  ## [6] "bestTune"     "call"         "dots"         "metric"       "control"
  ## [11] "finalModel"   "preProcess"   "trainingData" "resample"     "resampledCM"
  ## [16] "perfNames"    "maximize"     "yLimits"      "times"        "levels"

  expect_match(enFitMulti$method, "knn")
  ## enFitMulti$modelInfo
  expect_match(enFitMulti$modelType, "Classification")
  expect_match(typeof(enFitMulti$results), "list")
  expect_equal(enFitMulti$pred, NULL)
  ## enFitMulti$bestTune
  expect_equal(enFitMulti$bestTune$k, 7)
  ## enFitMulti$call
  expect_equal(enFitMulti$dots, list())
  expect_match(enFitMulti$metric, "Accuracy")
  ## enFitMulti$control
  ##enFitMulti$finalModel
  expect_equal(enFitMulti$preProcess, NULL)
  expect_equal(dim(enFitMulti$trainingData)[1], 1000)
  expect_equal(dim(enFitMulti$trainingData)[2], 3)
  expect_equal(dim(enFitMulti$resample)[1], 3)
  ## expect_equal(dim(enFitMulti$resample)[2], 16) 15
  expect_equal(dim(enFitMulti$resampledCM)[1], 6)
  expect_equal(dim(enFitMulti$resampledCM)[2], 11)
  ## expect_equal(length(enFitMulti$perfNames), 15)  ## 14
  expect_true(enFitMulti$maximize)
  expect_equal(enFitMulti$yLimits, NULL)
  ##enFitMulti$times
  ## enFitMulti$levels

})



### test 'feature_categories.R' functions used in 'generate_feature_table.R'

test_that("feature_categories.R functions", {

  maf = readRDS(system.file("extdata", "data_mutations_uniprot.rds", package = "MolecularDiagnosis"))
  seg = readRDS(system.file("extdata", "msk_impact_2017_data_cna_hg19_seg.rds", package = "MolecularDiagnosis"))
  cn = readRDS(system.file("extdata", "data_CNA.rds", package = "MolecularDiagnosis"))
  SV = readRDS(system.file("extdata", "data_fusions.rds", package = "MolecularDiagnosis"))
  ## mutations()
  ## mutations(maf = maf)
  expect_equal(dim(mutations(maf = head(maf, 500)))[1], 91)
  expect_equal(dim(mutations(maf = head(maf, 500)))[2], 343)
  ## truncating_mutations()
  ## truncating_mutations(maf = maf)
  expect_equal(dim(truncating_mutations(maf = head(maf, 500)))[1], 57)
  expect_equal(dim(truncating_mutations(maf = head(maf, 500)))[2], 342)
  ## hotspots()
  ## hotspots(maf = maf)
  expect_equal(dim(hotspots(maf = head(maf, 500)))[1], 69)
  expect_equal(dim(hotspots(maf = head(maf, 500)))[2], 4141)
  ## focal_cn_portal()
  ## focal_cn_portal(cn = cn)))
  expect_equal(dim(focal_cn_portal(cn = head(cn)))[1], 10945)
  expect_equal(dim(focal_cn_portal(cn = head(cn)))[2], 11)
  ## broad_cn()
  ## broad_cn(seg = seg)))
  expect_equal(dim(broad_cn(seg = head(seg)))[1], 1)
  expect_equal(dim(broad_cn(seg = head(seg)))[2], 85)
  ## fusions()
  ## fusions(SV = SV)))
  expect_equal(fusions(SV=head(SV))$EML4_ALK_fusion, 1)
  ## cn_burden()
  ## cn_burden(seg = seg)))
  expect_equal(cn_burden(seg = head(seg))$CN_Burden, 0.3)

})



maf = readRDS(system.file("extdata", "data_mutations_uniprot.rds", package = "MolecularDiagnosis"))

seg = readRDS(system.file("extdata", "msk_impact_2017_data_cna_hg19_seg.rds", package = "MolecularDiagnosis"))

cn = readRDS(system.file("extdata", "data_CNA.rds", package = "MolecularDiagnosis"))

SV = readRDS(system.file("extdata", "data_fusions.rds", package = "MolecularDiagnosis"))

## id_diagnosis = function(clinical = clinical, selected_cancer_types = "F", ps = "F", input_cancertypes = NULL){

## mutations = function(maf = maf){

## truncating_mutations = function(maf = maf) {

## hotspots <- function( maf = maf ){

## mutations = function(maf = maf){

## truncating_mutations = function(maf = maf) {

## hotspots <- function( maf = maf ){

## focal_cn_portal = function( cna = cn ) {

## focal_cn_AHD <- function( seg = seg ) {

## broad_cn = function(seg = seg, log_ratio_threshold = 0.2){

## fusions = function(SV = SV){

## structural_variants = function(clinical = clinical, SV = SV) {

## mutational_signatures = function( maf = maf){

## mutation_count = function(maf = maf){

## cn_burden = function(seg = seg, log_ratio_threshold = 0.2){

## gender = function(clinical = clinical){

## msi = function(clinical = clinical){

## KnownSignature = function( clinical = clinical ) {









