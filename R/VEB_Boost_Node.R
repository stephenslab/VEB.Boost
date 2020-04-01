#' @import data.tree

### Node Object ###

VEBBoostNode <- R6::R6Class(
  "VEBBoostNode",
  public = list(
    operator = NULL, # either "+" or "*" for internal nodes, NULL for terminal nodes

    currentFit = NULL, # current fit for fitting function

    fitFunction = NULL, # function that takes in predictors X, response Y, variances sigma2, and returns the fit
    # the fit must have fields mu1 (first moment), mu2 (second moment), KL_div (KL divergence from q to g),
    # and returns a function, predFunction(X_new, method = c(1, 2)) that takes in new data X and returns our prediction for the given moment

    predFunction = NULL, # function to predict based on current fit

    constCheckFunction = NULL, # function to check if fit is constant

    family = "gaussian",

    AddChildVEB = function(name, check = c("check", "no-warn", "no-check"), ...) { # add VEB node as child
      child = VEBBoostNode$new(as.character(name), check, ...)
      return(invisible(self$AddChildNode(child)))
    },

    updateMoments = function() { # after updating the fit, pass changes to moments up to parent node
      if (!self$isLeaf) { # if not at a leaf, update moments
        children_mu1 = sapply(self$children, function(x) x$mu1)
        children_mu2 = sapply(self$children, function(x) x$mu2)
        if (self$operator == "+") {
          private$.mu1 = children_mu1[, 1] + children_mu1[, 2]
          private$.mu2 = children_mu2[, 1] + children_mu2[, 2] + (2 * children_mu1[, 1] * children_mu1[, 2])
        } else {
          private$.mu1 = children_mu1[, 1] * children_mu1[, 2])
          private$.mu2 = children_mu2[, 1] * children_mu2[, 2])
        }
      }

      return(invisible(self))
    },

    updateMomentsAll = function() { # after updating the fit, pass changes to moments up to internal nodes
      if (!self$isLeaf) { # if not at a leaf, update moments
        children_mu1 = sapply(self$children, function(x) x$mu1)
        children_mu2 = sapply(self$children, function(x) x$mu2)
        if (self$operator == "+") {
          private$.mu1 = rowsums(children_mu1)
          private$.mu2 = rowsums(children_mu2) + 2*rowprods(children_mu1)
        } else {
          private$.mu1 = rowprods(children_mu1)
          private$.mu2 = rowprods(children_mu2)
        }
      }
      if (self$isRoot) { # if at root, stop
        return(invisible(self))
      } else { # else, update parents moments
        return(invisible(self$parent$updateMomentsAll()))
      }
    },
    
    updateCurrentInputs = function(currentInputs) { # update what current response and variances are (currentInputs = list(current_Y, current_sigma2)
      if (self$isRoot) {
        return(currentInputs)
      }
      if (self$parent$operator == "+") {
        currentInputs$Y = currentInputs$Y - self$siblings[[1]]$mu1
        currentInputs$sigma2 = currentInputs$sigma2
      } else {
        currentInputs$Y = currentInputs$Y * (self$siblings[[1]]$mu1 / self$siblings[[1]]$mu2)
        currentInputs$sigma2 = currentInputs$sigma2 / self$siblings[[1]]$mu2
      }
      return(currentInputs)
    }, 

    updateFit = function(currentInputs = NULL) { # function to update currentFit
      if (self$isLeaf) {
        if (is.null(currentInputs)) { # when starting at mu_0
          currentInputs = list(Y = self$Y, sigma2 = self$sigma2)
        } else { # else, update inputs
          currentInputs = self$updateCurrentInputs(currentInputs)
        }
        self$currentFit = self$fitFunction(X = self$X, Y = currentInputs$Y, sigma2 = currentInputs$sigma2, init = self$currentFit)
      }
      if (!self$isRoot) {
        self$parent$updateMoments()
        # now, revert to inputs at parent
        if (self$parent$operator == "+") {
          currentInputs$Y = currentInputs$Y + self$siblings[[1]]$mu1
          currentInputs$sigma2 = currentInputs$sigma2
        } else {
          currentInputs$Y = currentInputs$Y / (self$siblings[[1]]$mu1 / self$siblings[[1]]$mu2)
          currentInputs$sigma2 = currentInputs$sigma2 * self$siblings[[1]]$mu2
        }
      }
      return(currentInputs)
    },

    update_sigma2 = function() { # function to update sigma2
      self$sigma2 = ((sum(self$root$Y^2) - 2*sum(self$root$Y*self$root$mu1) + sum(self$root$mu2))) / length(self$root$Y)
      return(invisible(self))
    },

    convergeFit = function(tol = 1e-3, update_sigma2 = FALSE, update_ELBO_progress = TRUE, verbose = TRUE) {
      ELBOs = numeric(1000)
      ELBOs[1] = -Inf
      ELBOs[2] = self$root$ELBO
      i = 2
      while (abs(ELBOs[i] - ELBOs[i-1]) > tol) {
        currentInputs = NULL
        self$root$Do(function(x) currentInputs = x$updateFit(currentInputs), traversal = 'post-order')
        if (update_sigma2) {
          self$root$update_sigma2()
        }
        i = i+1
        if (i > length(ELBOs)) { # double size of ELBOs for efficiency rather than making it bigger each iteration
          ELBOs = c(ELBOs, rep(0, i))
        }
        ELBOs[i] = self$root$ELBO
        if (verbose & ((i %% 100) == 0)) {
          cat(paste("ELBO: ", ELBOs[i], sep = ""))
          cat("\n")
        }
      }
      ELBOs = ELBOs[2:i]

      if (update_ELBO_progress) {
        self$ELBO_progress = ELBOs
      }
      return(invisible(self$root))
    },

    AddSiblingVEB = function(learner, operator = c("+", "*"), combine_name) { # function to add subtree as sibling to given node, combining with operator
      # learner is tree to add to self
      # operator is how to combine them
      # name is what to call the combining node

      learner_copy = self$clone()
      learner_copy$X = NULL

      self$fitFunction = NULL
      self$currentFit = NULL
      self$predFunction = NULL
      self$constCheckFunction = NULL
      self$operator = operator
      self$name = combine_name

      self$AddChildNode(learner_copy)
      self$AddChildNode(learner)

      learner_copy$parent$updateMoments()

      return(invisible(self$root)) # return root, so we can assign entire tree to new value (in case we're splitting at root)

    },

    lockSelf = function(changeToConstant = TRUE) { # function node calls on itself to check if it's locked, and lock if necessary
      # if changeToConstant, change fitting function to constant
      if (self$isConstant) {
        self$isLocked = TRUE
        if (changeToConstant) {
          self$fitFunction = private$.fitFnConstComp
          self$predFunction = private$.predFnConstComp
          self$constCheckFunction = private$.constCheckFnConstComp
          self$updateFit()
          try({self$updateMomentsAll()}, silent = T) # needed to avoid "attempt to apply non-function" error
        }
      }
      return(invisible(self))
    },

    lockLearners = function(growMode = c("+", "*", "+*"), changeToConstant = TRUE) { # lock learners that should be locked
      # base_learners = Traverse(self$root, traversal = 'post-order', filterFun = function(x) x$isLeaf & !x$isLocked)
      # for (learner in base_learners) { # change any near-constant leaf to a constant and seal it off
      #   if (learner$isConstant) {
      #     learner$isLocked = TRUE
      #     learner$fitFunction = learner$.fitFnConstComp
      #     learner$predFunction = learner$.predFnConstComp
      #     learner$constCheckFunction = learner$.constCheckFnConstComp
      #     learner$updateFit()
      #     try({learner$updateMomentsAll()}, silent = T) # needed to avoid "attempt to apply non-function" error
      #   }
      # }
      self$root$Do(function(node) node$lockSelf(changeToConstant), filterFun = function(node) node$isLeaf & !node$isLocked, traversal ='post-order')

      base_learners = Traverse(self$root, filterFun = function(x) x$isLeaf & !x$isLocked)
      for (learner in base_learners) {
        if (learner$isRoot) { # if root, not locked, do this to avoid errors in next if statement
          next
        }
        if (learner$siblings[[1]]$isLocked && ((growMode %in% c("+", "*")) || learner$parent$siblings[[1]]$isLocked)) {
          learner$isLocked = TRUE
        }
      }

      return(invisible(self$root))
    },

    addLearnerAll = function(growMode = c("+", "*", "+*"), changeToConstant = TRUE) { # to each leaf, add a "+" and "*"
      self$root$lockLearners(growMode, changeToConstant) # lock learners

      base_learners = Traverse(self$root, filterFun = function(x) x$isLeaf & !x$isLocked)
      for (learner in base_learners) {
        fitFn = learner$fitFunction
        predFn = learner$predFunction
        constCheckFn = learner$constCheckFunction

        if (growMode %in% c("+", "*")) {
          learner_name = paste("mu_", learner$root$leafCount, sep = '')
          combine_name = paste("combine_", learner$root$leafCount, sep = '')

          add_fit = list(mu1 = rep(1 * (growMode == "*"), length(learner$Y)), mu2 = rep(1 * (growMode == "*"), length(learner$Y)), KL_div = 0)
          add_node = VEBBoostNode$new(learner_name, fitFunction = fitFn, predFunction = predFn, constCheckFunction = constCheckFn, currentFit = add_fit)
          learner$AddSiblingVEB(add_node, growMode, combine_name)
        } else {
          learner_name = paste("mu_", learner$root$leafCount, sep = '')
          combine_name = paste("combine_", learner$root$leafCount, sep = '')

          add_fit = list(mu1 = rep(0, length(learner$Y)), mu2 = rep(0, length(learner$Y)), KL_div = 0)
          add_node = VEBBoostNode$new(learner_name, fitFunction = fitFn, predFunction = predFn, constCheckFunction = constCheckFn, currentFit = add_fit)
          learner$AddSiblingVEB(add_node, "+", combine_name)

          learner_name = paste("mu_", learner$root$leafCount, sep = '')
          combine_name = paste("combine_", learner$root$leafCount, sep = '')

          mult_fit = list(mu1 = rep(1, length(learner$Y)), mu2 = rep(1, length(learner$Y)), KL_div = 0)
          mult_node = VEBBoostNode$new(learner_name, fitFunction = fitFn, predFunction = predFn, constCheckFunction = constCheckFn, currentFit = mult_fit)
          learner$children[[1]]$AddSiblingVEB(mult_node, "*", combine_name)
        }
      }

      return(invisible(self$root))
    },

    convergeFitAll = function(tol = 1e-3, update_sigma2 = FALSE, update_ELBO_progress = TRUE, growMode = "+*", changeToConstant = TRUE, verbose = FALSE) {
      self$addLearnerAll(growMode, changeToConstant)
      self$convergeFit(tol, update_sigma2, update_ELBO_progress, verbose)
      return(invisible(self$root))
    },

    predict.veb = function(X_new, moment = c(1, 2)) { # function to get prediction on new data
      self$root$Do(function(node) {
        if (1 %in% moment) {
          node$pred_mu1 = node$predFunction(X_new, node$currentFit, 1)
        }
        if (2 %in% moment) {
          node$pred_mu2 = node$predFunction(X_new, node$currentFit, 2)
        }
      }, filterFun = function(x) x$isLeaf)
      return(invisible(self$root))
    }

  ),
  private = list(
    .X = NULL, # predictors for node (if NA, use value of parent)
    .Y = NA, # response, only not NA at root
    .sigma2 = NA, # variance, only not NA at root
    .mu1 = NA, # current first moment
    .mu2 = NA, # current second moment
    .ELBO_progress = list(-Inf), # list of ELBOs for each iteration of growing and fitting the tree
    .pred_mu1 = NULL, # prediction based on predFunction and given new data (first moment)
    .pred_mu2 = NULL, # prediction based on predFunction and given new data (second moment)
    .isLocked = FALSE, # locked <=> V < V_tol, or both learners directly connected (sibling and parent's sibling) are constant
    .alpha = 0, # used in multi-class learner
    .ensemble = NULL, # used in multi-class learner, either reference to self, or multi-class learner object
    .fitFnConstComp = function(X, Y, sigma2, init) { # constant fit function
      if (length(sigma2) == 1) {
        sigma2 = rep(sigma2, length(Y))
      }
      intercept = weighted.mean(Y, 1/sigma2)
      KL_div = 0

      mu1 = intercept
      mu2 = intercept^2
      return(list(mu1 = mu1, mu2 = mu2, intercept = intercept, KL_div = KL_div))
    },
    .predFnConstComp = function(X_new, currentFit, moment = c(1, 2)) { # constant prediction function
      if (moment == 1) {
        return(currentFit$intercept)
      } else if (moment == 2) {
        return(currentFit$intercept^2)
      } else {
        stop("`moment` must be either 1 or 2")
      }
    },
    .constCheckFnConstComp = function(currentFit) {
      return(TRUE)
    },
    .g = function(x) { # inverse logit function
      1 / (1 + exp(-x))
    }
  ),

  active = list(

    ensemble = function(value) {
      if (missing(value)) {
        if (self$isRoot) {
          if (is.null(private$.ensemble)) {
            return(invisible(self))
          } else {
            return(private$.ensemble)
          }
        } else {
          return(invisible(self$root$ensemble))
        }
      }

      if (!self$isRoot) {
        stop("`$ensemble' cannot be modified except at the root", call. = FALSE)
      } else {
        private$.ensemble = value
      }
    },

    isEnsemble = function(value) { # TRUE if is ensemble (root AND not part of higher ensemble), else FALSE
      if (!missing(value)) {
        stop("`$isEnsemble' cannot be modified directly", call. = FALSE)
      }
      return((self$isRoot) && (is.null(private$.ensemble)))
    },

    mu1 = function(value) {
      if (!missing(value)) {
        stop("`$mu1` cannot be modified directly", call. = FALSE)
      }
      if (self$isLeaf) {
        mu1 = self$currentFit$mu1
      } else {
        mu1 = private$.mu1
      }
      if (length(mu1) == 1) {
        mu1 = rep(mu1, length(self$root$raw_Y))
      }
      return(mu1)
    },

    mu2 = function(value) {
      if (!missing(value)) {
        stop("`$mu2` cannot be modified directly", call. = FALSE)
      }
      if (self$isLeaf) {
        mu2 = self$currentFit$mu2
      } else {
        mu2 = private$.mu2
      }
      if (length(mu2) == 1) {
        mu2 = rep(mu2, length(self$root$raw_Y))
      }
      return(mu2)
    },

    KL_div = function(value) { # KL divergence from q to g of learner defined by sub-tree with this node as the root
      if (!missing(value)) {
        stop("`$KL_div` cannot be modified directly", call. = FALSE)
      }

      if (self$isLeaf) {
        return(self$currentFit$KL_div)
      }
      return(sum(sapply(self$children, function(x) x$KL_div)))
    },

    X = function(value) { # predictor for node
      if (missing(value)) {
        if (self$isRoot) {
          if (self$isEnsemble) {
            return(private$.X)
          } else {
            return(self$ensemble$X)
          }
        }
        if (is.null(private$.X)) {
          return(self$parent$X)
        } else {
          return(private$.X)
        }
      }
      private$.X = value
    },

    Y = function(value) { # response for sub-tree
      if (missing(value)) {
        if (self$isRoot) {
          if (self$root$family == "gaussian") {
            return(private$.Y)
          } else if (self$root$family == "binomial") {
            d = self$d
            return((private$.Y - .5 + (self$alpha * d)) / d) # alpha*d needed for multi-class case
          } else {
            stop("family must be either 'gaussian' or 'binomial")
          }
        }
        if (self$parent$operator == "+") {
          return(self$parent$Y - self$siblings[[1]]$mu1)
        }
        if (self$parent$operator == "*") {
          return(self$parent$Y * (self$siblings[[1]]$mu1 / self$siblings[[1]]$mu2))
        }
      } else {
        if (self$isRoot) {
          private$.Y = value
        } else {
          stop("`$Y` cannot be modified directly except at the root node", call. = FALSE)
        }
      }
    },

    raw_Y = function(value) { # raw value of private$.Y, only needed for logistic ELBO
      if (!missing(value)) {
        stop("`$raw_Y` cannot be modified directly", call. = FALSE)
      }
      if (self$isRoot) {
        return(private$.Y)
      } else {
        return(self$parent$raw_Y)
      }
    },

    sigma2 = function(value) { # variance for sub-tree
      if (missing(value)) {
        if (self$isRoot) {
          if (self$root$family == "gaussian") {
            return(private$.sigma2)
          } else if (self$root$family == "binomial") {
            return(1 / self$d)
          } else {
            stop("family must be either 'gaussian' or 'binomial")
          }
        }
        if (self$parent$operator == "+") {
          s2 = self$parent$sigma2
          if (length(s2) == 1) {
            s2 = rep(s2, length(self$raw_Y))
          }
          return(s2)
        }
        if (self$parent$operator == "*") {
          s2 = self$parent$sigma2 / self$siblings[[1]]$mu2
          if (length(s2) == 1) {
            s2 = rep(s2, length(self$raw_Y))
          }
          return(s2)
        }
      } else {
        if (self$isRoot) {
          private$.sigma2 = value
        } else {
          stop("`$sigma2` cannot be modified directly except at the root node", call. = FALSE)
        }
      }
    },

    ELBO_progress = function(value) { # ELBO progress of tree
      if (missing(value)) {
        if (self$isRoot) {
          return(private$.ELBO_progress)
        } else {
          return(self$root$ELBO_progress)
        }
      } else {
        if (self$isRoot) { # append ELBOs in list
          private$.ELBO_progress[[length(private$.ELBO_progress) + 1]] = value
        } else {
          stop("`$ELBO_progress` cannot be modified directly except at the root node", call. = FALSE)
        }
      }
    },

    ELBO = function(value) { # ELBO for entire tree
      if (!missing(value)) {
        stop("`$ELBO` cannot be modified directly", call. = FALSE)
      }
      # make s2 vector of variances
      #s2 = ifelse(length(self$sigma2) == 1, rep(self$sigma2, length(self$Y)), self$sigma2) # something weird w/ this if-else, not sure why
      if (self$root$family == "gaussian") {
        s2 = self$sigma2
        if (length(s2) == 1) {
          s2 = rep(s2, length(self$raw_Y))
        }
        return(
          (-.5 * sum(log(2*pi*s2))) -
            (.5 * (sum(self$Y^2 / s2) - 2*sum(self$Y * self$mu1 / s2) + sum(self$mu2 / s2))) -
            self$KL_div
        )
      } else if (self$root$family == "binomial") {
        d = self$root$d
        xi = self$root$xi
        return(
          sum(log(private$.g(xi))) + sum((xi / 2)*(d*xi - 1)) + sum((self$root$raw_Y - .5) * self$root$mu1) -
            .5*sum(self$root$mu2 * d) - self$KL_div
        )
      } else {
        stop("family must be either 'gaussian' or 'binomial")
      }
    },

    alpha = function(value) { # used in multi-class learner
      if (!missing(value)) {
        stop("`$alpha' cannot be modified directly except at the ensemble level", call. = FALSE)
      }
      if (self$isEnsemble) { # if isEnsemble, e.g. if we're in linear or logistic case, return root$.alpha (SHOULD BE 0 IN THIS CASE)
        return(private$.alpha)
      }
      return(self$ensemble$alpha) # otherwise, return ensemble's alpha
    },

    xi = function(value) { # optimal variational parameters, set to +sqrt(mu2)
      if (!missing(value)) {
        stop("`$xi` cannot be modified directly", call. = FALSE)
      }
      return(sqrt(self$root$mu2 + self$alpha^2 - 2*self$alpha*self$root$mu1))
    },

    d = function(value) { # d == 1/xi * (g(xi) - .5), n x K matrix
      if (!missing(value)) {
        stop("'$d' cannot be modified directly", call. = FALSE)
      }
      xi = self$xi
      g_xi = private$.g(xi) # matrix of g(xi_i,k), pre-compute once
      d = ((g_xi - .5) / xi) # matrix of (g(xi_i,k) - .5) / xi_i,k, pre-compute once
      d[xi == 0] = .25 # case of 0/0 (e.g. x_i is all 0), use L'Hopital
      return(d)
    },

    isConstant = function(value) { # uses node's constant check function and current fit to see if it is essentially a constant function
      if (!missing(value)) {
        stop("'$isConstant' cannot be modified directly", call. = FALSE)
      }
      return(self$constCheckFunction(self$currentFit))
    },

    isLocked = function(value) { # locked <=> V < V_tol, or both learners directly connected (sibling and parent's sibling) are constant
      if (missing(value)) {
        if (self$isLeaf) {
          return(private$.isLocked)
        } else {
          return(all(sapply(self$children, function(x) x$isLocked)))
        }
      } else {
        if (self$isLeaf) {
          private$.isLocked = value
        } else {
          stop("`$isLocked` cannot be modified directly except at leaf nodes", call. = FALSE)
        }
      }
    },

    pred_mu1 = function(value) { # predicted first moment given new data
      if (!missing(value)) {
        if (self$isLeaf) {
          private$.pred_mu1 = value
          return(invisible(self))
        } else {
          stop("`$pred_mu1` cannot be modified directly except at leaf nodes by calling `$predict.veb`", call. = FALSE)
        }
      }
      if (self$isLeaf) {
        return(private$.pred_mu1)
      } else {
        children_pred_mu1 = sapply(self$children, function(x) x$pred_mu1)
        if (self$operator == "+") {
          return(children_pred_mu1[, 1] + children_pred_mu1[, 2])
        } else {
          return(children_pred_mu1[, 1] * children_pred_mu1[, 2])
        }
      }
    },

    pred_mu2 = function(value) { # predicted second moment given new data
      if (!missing(value)) {
        if (self$isLeaf) {
          private$.pred_mu2 = value
          return(invisible(self))
        } else {
          stop("`$pred_mu2` cannot be modified directly except at leaf nodes by calling `$predict.veb`", call. = FALSE)
        }
      }
      if (self$isLeaf) {
        return(private$.pred_mu2)
      } else if (self$operator == "+") {
        children_pred_mu1 = sapply(self$children, function(x) x$pred_mu1)
        children_pred_mu2 = sapply(self$children, function(x) x$pred_mu2)
        return(children_pred_mu2[, 1] + children_pred_mu2[, 2] + (2 * children_pred_mu1[, 1] * children_pred_mu1[, 2]))
      } else {
        children_pred_mu2 = sapply(self$children, function(x) x$pred_mu2)
        return(children_pred_mu2[, 1] * children_pred_mu2[, 2])
      }
    }
  ),
  inherit = data.tree::Node,
  lock_objects = FALSE
)
