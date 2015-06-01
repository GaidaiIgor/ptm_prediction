# ---- functions
library(rpart)
library(e1071)
library(rpart.plot)
library(randomForest)
library(lattice)
library(gbm)
library(ROCR)

error_fun_auc_wrap = function(ordering)
{
  error_fun_auc = function(actual, predicted) 
  {
    pred = prediction(predicted, actual, ordering)
    -performance(pred, "auc")@y.values[[1]]
  }
  error_fun_auc
}

predict_func = function(x, newdata)
{
  predict(x, newdata)[, 1]
}

error_fun_auc = function(actual, predicted) 
{
  pred = prediction(predicted, actual, c("unmodified", "modified"))
  result = -performance(pred, "auc")@y.values[[1]]
  if (result > -0.5)
  {
    result = -1 - result
    #write("auc error", stderr())
  }
  result
}

plot_roc = function(predicted, actual, ordering = levels(actual))
{
  pred = prediction(predicted, actual, ordering)
  result = performance(pred, "auc")@y.values
  if (result < 0.5)
  {
    pred = prediction(predicted, actual)
    result = performance(pred, "auc")@y.values
  }
  plot(performance(pred, "tpr", "fpr"))
  result
}

get_random_partition = function(data, dependence_class, equally_distribute = 1 / test_share, test_share = 1 / 3)
{
  cv_folds = dependent_observation_cv_folds(dependence_class, 1 / test_share, equally_distribute)
  cv_folds = sample(cv_folds)
  test_ind = unlist(cv_folds[1])
  train_ind = unlist(cv_folds[-1])
  list(data[-test_ind, ], data[test_ind, ])
}

apply_to_columns = function(data, columns, func)
{
  data[, columns] = data.frame(lapply(data[, columns], func))
  return (data)
}

prepare_data = function(data, select_columns = NULL)
{
  data = na.omit(data)
  factor_prefix = c("sec_struct")
  abs_prefix = c("omega")
  # remove trailing numbers
  unified_names = gsub("\\d+$", "", names(data))
  columns_to_process =  unified_names %in% factor_prefix
  data = apply_to_columns(data, columns_to_process, factor)
  columns_to_process = unified_names %in% abs_prefix
  data = apply_to_columns(data, columns_to_process, abs)
  if (!is.null(select_columns))
  {
    data = data[, unified_names %in% select_columns]
  }
  return (data)
}

expand_formula = function(formula, data)
{
  reformulate(attr(terms(formula, data = data), "term.labels"), formula[[2]])
}

new_tune = function (method, train.x, train.y = NULL, data = list(), validation.x = NULL, 
                    validation.y = NULL, ranges = NULL, predict.func = predict, 
                    tunecontrol = tune.control(), dependence_class = NULL, 
                    equally_distribute = tunecontrol$cross, ...) 
{
  call <- match.call()
  resp <- function(formula, data) {
    model.response(model.frame(formula, data))
  }
  classAgreement <- function(tab) {
    n <- sum(tab)
    if (!is.null(dimnames(tab))) {
      lev <- intersect(colnames(tab), rownames(tab))
      p0 <- sum(diag(tab[lev, lev]))/n
    }
    else {
      m <- min(dim(tab))
      p0 <- sum(diag(tab[1:m, 1:m]))/n
    }
    p0
  }
  if (tunecontrol$sampling == "cross") 
    validation.x <- validation.y <- NULL
  useFormula <- is.null(train.y)
  if (useFormula && (is.null(data) || length(data) == 0)) 
    data <- model.frame(train.x)
  if (is.vector(train.x)) 
    train.x <- t(t(train.x))
  if (is.data.frame(train.y)) 
    train.y <- as.matrix(train.y)
  if (!is.null(validation.x)) 
    tunecontrol$fix <- 1
  n <- nrow(if (useFormula) 
    data
    else train.x)
  perm.ind <- sample(n)
  if (tunecontrol$sampling == "cross") {
    if (tunecontrol$cross > n) 
      stop(sQuote("cross"), " must not exceed sampling size!")
    if (tunecontrol$cross == 1) 
      stop(sQuote("cross"), " must be greater than 1!")
  }
  train_set = if (useFormula) data else train.x
  train.ind <- if (tunecontrol$sampling == "cross") 
    if (is.null(dependence_class))
      tapply(1:n, cut(1:n, breaks = tunecontrol$cross), function(x) perm.ind[-x])
    else 
    {
      cv_folds = dependent_observation_cv_folds(dependence_class, tunecontrol$cross, equally_distribute)
      cv_folds_train = Map(function(x) seq_len(length(dependence_class))[-x], cv_folds)
      cv_folds_train
    }
      
  else if (tunecontrol$sampling == "fix") 
    list(perm.ind[1:trunc(n * tunecontrol$fix)])
  else lapply(1:tunecontrol$nboot, function(x) sample(n, n * 
                                                        tunecontrol$boot.size, replace = TRUE))
  parameters <- if (is.null(ranges)) 
    data.frame(dummyparameter = 0)
  else expand.grid(ranges)
  p <- nrow(parameters)
  if (!is.logical(tunecontrol$random)) {
    if (tunecontrol$random < 1) 
      stop("random must be a strictly positive integer")
    if (tunecontrol$random > p) 
      tunecontrol$random <- p
    parameters <- parameters[sample(1:p, tunecontrol$random), 
                             ]
  }
  model.variances <- model.errors <- c()
  for (para.set in 1:p) {
    sampling.errors <- c()
    for (sample in 1:length(train.ind)) {
      repeat.errors <- c()
      for (reps in 1:tunecontrol$nrepeat) {
        pars <- if (is.null(ranges)) 
          NULL
        else lapply(parameters[para.set, , drop = FALSE], 
                    unlist)
        model <- if (useFormula) 
          do.call(method, c(list(train.x, data = data, 
                                 subset = train.ind[[sample]]), pars, list(...)))
        else do.call(method, c(list(train.x[train.ind[[sample]], 
                                            ], y = train.y[train.ind[[sample]]]), pars, 
                               list(...)))
        pred <- predict.func(model, if (!is.null(validation.x)) 
          validation.x
          else if (useFormula) 
            data[-train.ind[[sample]], , drop = FALSE]
          else if (inherits(train.x, "matrix.csr")) 
            train.x[-train.ind[[sample]], ]
          else train.x[-train.ind[[sample]], , drop = FALSE])
        true.y <- if (!is.null(validation.y)) 
          validation.y
        else if (useFormula) {
          if (!is.null(validation.x)) 
            resp(train.x, validation.x)
          else resp(train.x, data[-train.ind[[sample]], 
                                  ])
        }
        else train.y[-train.ind[[sample]]]
        if (is.null(true.y)) 
          true.y <- rep(TRUE, length(pred))
        repeat.errors[reps] <- if (!is.null(tunecontrol$error.fun)) 
          tunecontrol$error.fun(true.y, pred)
        else if ((is.logical(true.y) || is.factor(true.y)) && 
                   (is.logical(pred) || is.factor(pred) || is.character(pred))) 
          1 - classAgreement(table(pred, true.y))
        else if (is.numeric(true.y) && is.numeric(pred)) 
          crossprod(pred - true.y)/length(pred)
        else stop("Dependent variable has wrong type!")
      }
      sampling.errors[sample] <- tunecontrol$repeat.aggregate(repeat.errors)
    }
    model.errors[para.set] <- tunecontrol$sampling.aggregate(sampling.errors)
    model.variances[para.set] <- tunecontrol$sampling.dispersion(sampling.errors)
  }
  best <- which.min(model.errors)
  pars <- if (is.null(ranges)) 
    NULL
  else lapply(parameters[best, , drop = FALSE], unlist)
  structure(list(best.parameters = parameters[best, , drop = FALSE], 
                 best.performance = model.errors[best], method = if (!is.character(method)) deparse(substitute(method)) else method, 
                 nparcomb = nrow(parameters), train.ind = train.ind, sampling = switch(tunecontrol$sampling, 
                                                                                       fix = "fixed training/validation set", bootstrap = "bootstrapping", 
                                                                                       cross = if (tunecontrol$cross == n) "leave-one-out" else paste(tunecontrol$cross, 
                                                                                                                                                      "-fold cross validation", sep = "")), performances = if (tunecontrol$performances) cbind(parameters, 
                                                                                                                                                                                                                                               error = model.errors, dispersion = model.variances), 
                 best.model = if (tunecontrol$best.model) {
                   modeltmp <- if (useFormula) do.call(method, c(list(train.x, 
                                                                      data = data), pars, list(...))) else do.call(method, 
                                                                                                                   c(list(x = train.x, y = train.y), pars, list(...)))
                   call[[1]] <- as.symbol("best.tune")
                   modeltmp$call <- call
                   modeltmp
                 }), class = "tune")
}

gbm_wrap = function(..., data, subset = seq_len(nrow(data)), dependence_gbm) 
{
  new_gbm(..., data = data[subset, ], dependence_gbm = dependence_gbm[subset])
}

# my.predict.gbm <- function(x, newdata, ...) {
#   pred <- predict(x, newdata)
#   labels <- dimnames(pred)[[2]]
#   dim(pred) <- dim(pred)[1:2]
#   pred <- factor(max.col(pred), levels = seq_along(labels),
#                  labels = labels)
#   pred
# }

# my_predict_gbm = function(model, new_data) 
# {
#   pred = predict(model, new_data, type = "response")[, , 1]
#   colnames(pred)[max.col(pred)]
# }

new_gbm = function (formula = formula(data), distribution = "bernoulli", 
          data = list(), weights, var.monotone = NULL, n.trees = 100, 
          interaction.depth = 1, n.minobsinnode = 10, shrinkage = 0.001, 
          bag.fraction = 0.5, train.fraction = 1, cv.folds = 0, keep.data = TRUE, 
          verbose = "CV", class.stratify.cv = NULL, n.cores = NULL, 
          dependence_gbm = NULL, equally_distribute = cv.folds) 
{
  theCall <- match.call()
  lVerbose <- if (!is.logical(verbose)) {
    FALSE
  }
  else {
    verbose
  }
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "weights", "offset"), names(mf), 
             0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.pass
  mf[[1]] <- as.name("model.frame")
  m <- mf
  mf <- eval(mf, parent.frame())
  Terms <- attr(mf, "terms")
  y <- model.response(mf)
  if (missing(distribution)) {
    distribution <- guessDist(y)
  }
  else if (is.character(distribution)) {
    distribution <- list(name = distribution)
  }
  w <- model.weights(mf)
  offset <- model.offset(mf)
  var.names <- attributes(Terms)$term.labels
  x <- model.frame(terms(reformulate(var.names)), data, na.action = na.pass)
  response.name <- as.character(formula[[2]])
  lVerbose <- if (!is.logical(verbose)) {
    FALSE
  }
  else {
    verbose
  }
  class.stratify.cv <- getStratify(class.stratify.cv, distribution)
  group <- NULL
  num.groups <- 0
  if (distribution$name != "pairwise") {
    nTrain <- floor(train.fraction * nrow(x))
  }
  else {
    distribution.group <- distribution[["group"]]
    if (is.null(distribution.group)) {
      stop("For pairwise regression, the distribution parameter must be a list with a parameter 'group' for the a list of the column names indicating groups, for example list(name=\"pairwise\",group=c(\"date\",\"session\",\"category\",\"keywords\")).")
    }
    i <- match(distribution.group, colnames(data))
    if (any(is.na(i))) {
      stop("Group column does not occur in data: ", distribution.group[is.na(i)])
    }
    group <- factor(do.call(paste, c(data[, distribution.group, 
                                          drop = FALSE], sep = ":")))
    if ((!missing(weights)) && (!is.null(weights))) {
      w.min <- tapply(w, INDEX = group, FUN = min)
      w.max <- tapply(w, INDEX = group, FUN = max)
      if (any(w.min != w.max)) {
        stop("For distribution 'pairwise', all instances for the same group must have the same weight")
      }
      w <- w * length(w.min)/sum(w.min)
    }
    perm.levels <- levels(group)[sample(1:nlevels(group))]
    group <- factor(group, levels = perm.levels)
    ord.group <- order(group, -y)
    group <- group[ord.group]
    y <- y[ord.group]
    x <- x[ord.group, , drop = FALSE]
    w <- w[ord.group]
    num.groups.train <- max(1, round(train.fraction * nlevels(group)))
    nTrain <- max(which(group == levels(group)[num.groups.train]))
    Misc <- group
  }
  cv.error <- NULL
  if (cv.folds > 1) {
    cv.results <- new_gbm_cross_val(cv.folds, nTrain, n.cores, 
                              class.stratify.cv, data, x, y, offset, distribution, 
                              w, var.monotone, n.trees, interaction.depth, n.minobsinnode, 
                              shrinkage, bag.fraction, var.names, response.name, 
                              group, dependence_gbm, equally_distribute)
    cv.error <- cv.results$error
    p <- cv.results$predictions
  }
  gbm.obj <- gbm.fit(x, y, offset = offset, distribution = distribution, 
                     w = w, var.monotone = var.monotone, n.trees = n.trees, 
                     interaction.depth = interaction.depth, n.minobsinnode = n.minobsinnode, 
                     shrinkage = shrinkage, bag.fraction = bag.fraction, nTrain = nTrain, 
                     keep.data = keep.data, verbose = lVerbose, var.names = var.names, 
                     response.name = response.name, group = group)
  gbm.obj$train.fraction <- train.fraction
  gbm.obj$Terms <- Terms
  gbm.obj$cv.error <- cv.error
  gbm.obj$cv.folds <- cv.folds
  gbm.obj$call <- theCall
  gbm.obj$m <- m
  if (cv.folds > 0) {
    gbm.obj$cv.fitted <- p
  }
  if (distribution$name == "pairwise") {
    gbm.obj$ord.group <- ord.group
    gbm.obj$fit <- gbm.obj$fit[order(ord.group)]
  }
  return(gbm.obj)
}

new_gbm_cross_val = function(cv.folds, nTrain, n.cores, class.stratify.cv, data, 
          x, y, offset, distribution, w, var.monotone, n.trees, interaction.depth, 
          n.minobsinnode, shrinkage, bag.fraction, var.names, response.name, 
          group, dependence_class = NULL, equally_distribute = cv.folds) 
{
  i.train <- 1:nTrain
  cv.group <- new_get_cv_group(distribution, class.stratify.cv, y, 
                         i.train, cv.folds, group, dependence_class, equally_distribute)
#   cv.group <- getCVgroup(distribution, class.stratify.cv, y, 
#                                i.train, cv.folds, group)
  cv.models <- gbmCrossValModelBuild(cv.folds, cv.group, n.cores, 
                                     i.train, x, y, offset, distribution, w, var.monotone, 
                                     n.trees, interaction.depth, n.minobsinnode, shrinkage, 
                                     bag.fraction, var.names, response.name, group)
  cv.error <- gbmCrossValErr(cv.models, cv.folds, cv.group, 
                             nTrain, n.trees)
  best.iter.cv <- which.min(cv.error)
  predictions <- gbmCrossValPredictions(cv.models, cv.folds, 
                                        cv.group, best.iter.cv, distribution, data[i.train, ], 
                                        y)
  list(error = cv.error, predictions = predictions)
}

flat_cv_folds = function(cv_folds, vector_size)
{
  result = vector("numeric", vector_size)
  for (i in seq_len(length(cv_folds)))
  {
    cv_fold = cv_folds[[i]]
    for (index in cv_fold)
    {
      result[index] = i
    }
  }
  
  result
}

# distributes dependent groups of observations among cv folds
# returns list of vectors with indexes of observations which belong to that fold
dependent_observation_cv_folds = function(dependence_class, cv_folds_number, 
                                              equally_distribute = cv_folds_number)
{
  classes = unique(dependence_class)
  dependence_class = data.frame(class = dependence_class, index = seq_len(length(dependence_class)))
  observation_groups = aggregate(index ~ class, dependence_class, c, simplify = FALSE)
  groups_size = unlist(Map(length, observation_groups[, 2]))
  observation_groups = observation_groups[order(-groups_size), ]
  heavy_groups = head(observation_groups, equally_distribute)
  light_groups = tail(observation_groups, nrow(observation_groups) - equally_distribute)
  
  cv_folds = replicate(cv_folds_number, vector("numeric"))
  for (observation_group in heavy_groups[, 2])
  {
    min_index = which.min(unlist(Map(length, cv_folds)))
    cv_folds[[min_index]] = c(cv_folds[[min_index]], observation_group)
  }
  fold_numbers = sample(rep(seq_len(cv_folds_number), length = nrow(light_groups)))
  for (i in seq_len(nrow(light_groups)))
  {
    cv_folds[[fold_numbers[i]]] = c(cv_folds[[fold_numbers[i]]], unlist(light_groups[i, 2]))
  }
  cv_folds
}

new_get_cv_group = function (distribution, class.stratify.cv, y, i.train, cv.folds, 
                          group, dependence_class = NULL, equally_distribute = cv.folds) 
{
  if (!is.null(dependence_class))
  {
    cv_folds = dependent_observation_cv_folds(dependence_class, cv.folds, equally_distribute)
    cv.group = flat_cv_folds(cv_folds, length(dependence_class))
  }
  else if (distribution$name %in% c("bernoulli", "multinomial") & 
        class.stratify.cv) {
    nc <- table(y[i.train])
    uc <- names(nc)
    if (min(nc) < cv.folds) {
      stop(paste("The smallest class has only", min(nc), 
                 "objects in the training set. Can't do", cv.folds, 
                 "fold cross-validation."))
    }
    cv.group <- vector(length = length(i.train))
    for (i in 1:length(uc)) {
      cv.group[y[i.train] == uc[i]] <- sample(rep(1:cv.folds, 
                                                  length = nc[i]))
    }
  }
  else if (distribution$name == "pairwise") {
    s <- sample(rep(1:cv.folds, length = nlevels(group)))
    cv.group <- s[as.integer(group[i.train])]
  }
  else {
    cv.group <- sample(rep(1:cv.folds, length = length(i.train)))
  }
  cv.group
}