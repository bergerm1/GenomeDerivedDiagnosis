get_fc <- function(rfe_model, feature_table, mClass = NULL){
  rf_object <- rfe_model$fit$finalModel$classifier
  dataT <- as.data.frame(rfe_model$fit$trainingData)
  rownames(dataT) <- names(rf_object$predicted)

  #Calculate feature contributions
  li <- getLocalIncrements(rf_object, dataT)
  log <- capture.output({
    fc <- featureContributions(rf_object, li, feature_table, mClass = mClass)
  }, file = "/dev/null")

  rownames(fc) <- feature_table$SAMPLE_ID
  as.data.table(fc, keep.rownames = "SAMPLE_ID")
}

predict_table <- function(
  rfe_model,
  feature_table,
  model_function_list,
  N_top_cancer_types = 3
){

  Cancer_Types <- rfe_model$obsLevels
  vars <- rfe_model$optVariables

  x <- feature_table[, vars, with = F]
  rownames(x) <- feature_table$SAMPLE_ID

  pred_mat <- as.data.table(matrix(nrow = nrow(x), model_function_list$prob(modelFit = rfe_model$fit$finalModel, newdata = x)))
  names(pred_mat) <- Cancer_Types
  pred <- cbind(feature_table[, .(SAMPLE_ID, Classification_Category, Cancer_Type)],
                pred_mat)
  pred[, paste0("Pred", 1:N_top_cancer_types) := as.list(
    Cancer_Types[head(n=N_top_cancer_types,
                      order(.SD,
                            decreasing = T))]),
    SAMPLE_ID,
    .SDcols = Cancer_Types]
  pred[, paste0("Conf", 1:N_top_cancer_types) := as.list(round(
    digits = 4,
    as.numeric(
      head(n=N_top_cancer_types,
           sort(.SD, decreasing = T))))),
    SAMPLE_ID,
    .SDcols = Cancer_Types]

  setcolorder(
    pred,
    c(
      "SAMPLE_ID",
      "Cancer_Type",
      "Classification_Category",
      paste0(c("Pred", "Conf"),
             rep(1:N_top_cancer_types, each = 2)),
      Cancer_Types
    )
  )
  pred_colnames <- c(
    "SAMPLE_ID",
    "Cancer_Type",
    "Classification_Category",
    paste0(c("Pred", "Conf"),
           rep(1:N_top_cancer_types, each = 2)),
    Cancer_Types
  )
  pred <- pred[, pred_colnames, with = F]

  pred
}

get_xval_predictions <- function(rfe_model, N_top_cancer_types = 3){
  Cancer_Types <- rfe_model$obsLevels

  pred_train <- setDT(rfe_model$fit$pred)
  setkeyv(pred_train, names(rfe_model$fit$bestTune))
  pred_train <- pred_train[rfe_model$fit$bestTune]
  pred_train <- pred_train[order(rowIndex)]

  pred_train[, paste0("Pred", 1:N_top_cancer_types) := as.list(
    Cancer_Types[head(n=N_top_cancer_types,
                      order(.SD,
                            decreasing = T))]),
    rowIndex,
    .SDcols = Cancer_Types]
  pred_train[, paste0("Conf", 1:N_top_cancer_types) := as.list(round(
    digits = 4,
    as.numeric(
      head(n=N_top_cancer_types,
           sort(.SD, decreasing = T))))),
    rowIndex,
    .SDcols = Cancer_Types]
  pred_train
}

.predict_new <- function(
  rfe_model,
  feature_table,
  model_function_list,
  N_top_cancer_types = 3,
  all_cancer_types = TRUE,
  N_top_variables = 10,
  all_variables = FALSE
){

  feature_table_lab <- feature_table[Classification_Category == "train", 1:8, with = F]

  Cancer_Types <- rfe_model$obsLevels
  vars <- rfe_model$optVariables

  x <- feature_table[order(SAMPLE_ID)][, vars, with = F]
  rownames(x) <- feature_table[order(SAMPLE_ID)]$SAMPLE_ID

  pred_train <- get_xval_predictions(rfe_model, N_top_cancer_types = 3)

  pred_train <- cbind(
    feature_table[Classification_Category == "train", 1:8, with = F],
    pred_train)
  pred_colnames <- c(
    "SAMPLE_ID",
    "Cancer_Type",
    "Classification_Category",
    paste0(c("Pred", "Conf"),
           rep(1:N_top_cancer_types, each = 2)),
    Cancer_Types
  )
  pred_train <- pred_train[, pred_colnames, with = F]

  pred_new <- predict_table(
    rfe_model,
    feature_table[Classification_Category != "train"],
    model_function_list = model_function_list,
    N_top_cancer_types = N_top_cancer_types
  )

  pred <- rbind(pred_train, pred_new)

  feature_table_pred <- merge(
    feature_table,
    pred[, .(SAMPLE_ID, Pred1)],
    by = "SAMPLE_ID"
  )

  fc <- feature_table_pred[, get_fc(rfe_model, .SD, mClass = Pred1), Pred1]

  fc <- merge(fc, feature_table[, .(SAMPLE_ID, Classification_Category)], by = "SAMPLE_ID")
  fc <- fc[order(SAMPLE_ID)]

  feature_present <- copy(x)
  for(var in c("Gender_F", "LogSNV_Mb", "LogINDEL_Mb", "CN_Burden")){
    if(var %in% names(feature_present))
      feature_present[, (var) := 1]
  }

  # all(feature_present == 0 | feature_present == 1)

  fc_sum <- rowSums(fc[, vars, with = F])

  scale_factor <- pred$Conf1 / fc_sum
  fc[, (vars) := as.data.table(.SD * scale_factor), .SDcols = vars]

  dummy_val <- 1000

  fc_present <- fc[, .SD + feature_present * dummy_val, .SDcols = vars]  ## label features that are present by adding 1000
  fc_present <- cbind(
    "SAMPLE_ID" = fc[[1]],
    fc_present
  )

  var_vec <- as.numeric(fc_present[1, vars, with = F])
  names(var_vec) <- names(fc_present[1, vars, with = F])
  # varimp_names <-
  top_imp_val <- fc_present[, {
    var_vec <- as.numeric(.SD)
    names(var_vec) <- names(.SD)
    names(var_vec) <- ifelse(var_vec > dummy_val / 2,
                             names(var_vec),
                             paste0("Absence of ", names(var_vec)))
    var_vec <- ifelse(var_vec > dummy_val / 2,
                      var_vec - dummy_val,
                      var_vec)
    var_vec_top <- head(n=N_top_variables,
                        sort(var_vec[!is.infinite(var_vec)], decreasing = T))
    out_list <- c(
      as.list(
        c(
          names(var_vec_top),
          round(digits = 6,
                as.numeric(var_vec_top))
        ))
      # ,list(list(as.list(var_vec)))
    )
    names(out_list) <- c(paste0("Var", 1:N_top_variables), paste0("Imp", 1:N_top_variables)
                         #, "VarList"
    )
    out_list
  },
  SAMPLE_ID,
  .SDcols = vars]
  top_imp_val_cols <- paste0("Imp", 1:N_top_variables)
  top_imp_val[, (top_imp_val_cols) := lapply(.SD, as.numeric), .SDcols = top_imp_val_cols]

  # fc[, ImpSum := fc_sum]

  setcolorder(
    top_imp_val,
    c(
      "SAMPLE_ID",
      # "ImpSum",
      paste0(c("Var", "Imp"),
             rep(1:N_top_variables, each = 2))
    )
  )
  pred_fc <- merge(pred, top_imp_val, by = "SAMPLE_ID")
  if(all_variables == TRUE) pred_fc <- merge(pred_fc, fc_present, by = "SAMPLE_ID")

  pred_fc
}

#' Predict new samples
#'
#' @param arg_line
#'
#' @return
#' @export
#' @import rfFC
predict_new <- function(arg_line = NA){

  ### process args
  if(!is.na(arg_line) | interactive()) {
    # print("reading from arg_line")
    raw_args <- unlist(stringr::str_split(arg_line, " "))
  } else {
    # print("batch")
    raw_args <- commandArgs(TRUE)
  }
  # print(raw_args)

  option_list <- list(
    optparse::make_option(c("-t", "--test"),
                          default = FALSE,
                          action="store_true",
                          help="train a limited test classifier")
  )
  if(any(sapply(
    option_list,
    function(option){
      option@short_flag == "-g"
    }))){
    stop("cannot use short option '-g', conflicts with Rscript --gui")
  }

  opts <- optparse::parse_args(
    optparse::OptionParser(option_list=option_list),
    args = raw_args,
    positional_arguments = TRUE
  )

  model_filename <- opts$args[[1]]; opts$args <- opts$args[-1]
  feature_table_filename <- opts$args[[1]]; opts$args <- opts$args[-1]
  prediction_table_filename <- opts$args[[1]]; opts$args <- opts$args[-1]

  rfe_model <- readRDS(model_filename)
  feature_table <- fread(feature_table_filename)

  pred <- .predict_new(
    rfe_model = rfe_model,
    feature_table = feature_table,
    model_function_list = MolecularDiagnosis::rfcal,
    N_top_cancer_types = 3,
    all_cancer_types = TRUE,
    N_top_variables = 10,
    all_variables = FALSE
  )

  openxlsx::write.xlsx(
    pred,
    prediction_table_filename
  )
}
