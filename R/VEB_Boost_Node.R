#' @import data.tree
#' @import R6

### Node Object ###

VEBBoostNode <- R6Class(
  "VEBBoostNode",
  public = list(
    operator = NULL, # either "+" or "*" for internal nodes, NULL for terminal nodes

    learner = NULL, # list containing fitFunction, predFunction, constCheckFunction, currentFit, X, X_test, growMode, and changeToConstant

    family = "gaussian",

    cutpoints = NULL, # cutpoints for ordinal regression

    weights = 1, # observation weights (constrained to be >0, mean of 1)

    my_class_index = NULL, # for multinomial.titsias

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
          private$.mu1 = children_mu1[, 1] * children_mu1[, 2]
          private$.mu2 = children_mu2[, 1] * children_mu2[, 2]
        }
      }

      return(invisible(self))
    },

    updateMomentsAll = function() { # after updating the fit, pass changes to moments up to internal nodes
      self$updateMoments()
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
      if (currentInputs$name == self$parent$name) {
        if (self$parent$operator == "+") {
          currentInputs$Y = currentInputs$Y - self$siblings[[1]]$mu1
          currentInputs$sigma2 = currentInputs$sigma2
        } else {
          currentInputs$Y = currentInputs$Y * (self$siblings[[1]]$mu1 / self$siblings[[1]]$mu2)
          currentInputs$sigma2 = currentInputs$sigma2 / self$siblings[[1]]$mu2
        }
        currentInputs$name = self$name
        return(currentInputs)
      } else {
        return(self$updateCurrentInputs(self$parent$updateCurrentInputs(currentInputs)))
      }
    },

    updateFit = function(currentInputs = NULL) { # function to update currentFit
      if (self$isLeaf) {
        if (is.null(currentInputs)) { # when starting at mu_0
          currentInputs = list(Y = self$Y, sigma2 = self$sigma2, name = self$name)
        } else { # else, update inputs
          currentInputs = self$updateCurrentInputs(currentInputs)
        }
        if (!self$isLocked && any(is.infinite(currentInputs$sigma2))) { # if somehow variances become infinite, just a safeguard....
          self$isLocked = TRUE
          self$learner = constLearner
        }

        self$learner$currentFit = self$learner$fitFunction(self$X, currentInputs$Y, currentInputs$sigma2, self$learner$currentFit)
      }
      self$updateMoments()
      if (!self$isRoot) {
        currentInputs$name = self$parent$name
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

    updateSigma2 = function() { # function to update sigma2
      if (self$root$family == "aft.loglogistic") {
        a1 = -sum(self$root$weights[attr(self$root$raw_Y, 'not_cens')])
        a2 = .5*(sum(self$root$weights[attr(self$root$raw_Y, 'left_cens')] * (attr(self$root$raw_Y, 'log_Y')[attr(self$root$raw_Y, 'left_cens'), 2] - self$root$mu1[attr(self$root$raw_Y, 'left_cens')])) +
                   sum(self$root$weights[attr(self$root$raw_Y, 'right_cens')] * (self$root$mu1[attr(self$root$raw_Y, 'right_cens')] - attr(self$root$raw_Y, 'log_Y')[attr(self$root$raw_Y, 'right_cens'), 1])) -
                   sum(self$root$weights[attr(self$root$raw_Y, 'int_cens')] * (attr(self$root$raw_Y, 'log_Y')[attr(self$root$raw_Y, 'int_cens'), 2] - attr(self$root$raw_Y, 'log_Y')[attr(self$root$raw_Y, 'int_cens'), 1])))
        a3 = -.5*sum(cbind(self$root$weights, self$root$weights)*self$root$d*(cbind(self$root$mu2, self$root$mu2) - 2*cbind(self$root$mu1, self$root$mu1)*attr(self$root$raw_Y, 'log_Y') + attr(self$root$raw_Y, 'log_Y')^2), na.rm = TRUE)
        nll.logscale = function(lS) {
          -(a1*lS + a2*exp(-lS) + a3*exp(-2*lS) +
              sum(attr(self$root$raw_Y, 'log_Y')[attr(self$root$raw_Y, 'int_cens'), 2])/exp(lS) +
              sum(log(-expm1((attr(self$root$raw_Y, 'log_Y')[attr(self$root$raw_Y, 'int_cens'), 1] - attr(self$root$raw_Y, 'log_Y')[attr(self$root$raw_Y, 'int_cens'), 2])/exp(lS)))))
        }
        lS = optim(par = .5*log(self$raw_sigma2), fn = nll.logscale, method = 'Brent', lower = -15, upper = 35)$par
        self$sigma2 = exp(2*lS)
      } else if (self$root$family == "ordinal.logistic") {
        d = self$root$d * cbind(self$root$weights, self$root$weights)
        n_k = aggregate(self$root$weights, by = list(self$root$raw_Y), FUN = sum)[, 2] #n_k = table(self$root$raw_Y)
        sum_d_k_2 = .5*do.call(rbind, lapply(1:max(self$root$raw_Y), function(k) colSums(d[which(self$root$raw_Y == k), ])))
        sum_dT_k = do.call(rbind, lapply(1:max(self$root$raw_Y), function(k) colSums(sweep(d[which(self$root$raw_Y == k), ], 1, self$root$mu1[which(self$root$raw_Y == k)], "*"))))
        fn = function(theta) {
          res = 0
          for (k in 1:length(theta)) {
            res = res + (theta[k]^2)*(sum_d_k_2[k, 2] + sum_d_k_2[k+1, 1]) - theta[k]*(.5*(n_k[k] - n_k[k+1]) + sum_dT_k[k, 2] + sum_dT_k[k+1, 1])
            if (k != 1) {
              res = res - n_k[k]*log(1 - exp(theta[k-1] - theta[k]))
            }
          }
          return(res)
        }
        gr = function(theta) {
          res = numeric(length(theta))
          for (k in 1:length(theta)) {
            res[k] = 2*theta[k]*(sum_d_k_2[k, 2] + sum_d_k_2[k+1, 1]) - (.5*(n_k[k] - n_k[k+1]) + sum_dT_k[k, 2] + sum_dT_k[k+1, 1])
            if (k != 1) {
              res[k] = res[k] - n_k[k]*(exp(theta[k] - theta[k-1]) - 1)^(-1)
            }
            if (k != length(theta)) {
              res[k] = res[k] + n_k[k+1]*(exp(theta[k+1] - theta[k]) - 1)^(-1)
            }
          }
          return(res)
        }

        ui = toeplitz(c(-1, 1, rep(0, length(n_k)-3)))
        ui[lower.tri(ui, diag = FALSE)] = 0
        ui = ui[-nrow(ui), , drop = FALSE]

        init = self$root$cutpoints
        attributes(init) = NULL
        theta = constrOptim(theta = init, f = fn, grad = gr, ui = ui, ci = rep(0, nrow(ui)))$par

        self$cutpoints = theta
        attr(self$cutpoints, "log_Y") = matrix(NA, nrow = length(self$root$raw_Y), ncol = 2)
        attr(self$cutpoints, "log_Y")[attr(self$root$raw_Y, 'left_cens'), 2] = self$cutpoints[1]
        for (k in 2:(length(self$cutpoints))) {
          attr(self$cutpoints, "log_Y")[which(self$root$raw_Y == k), 1] = self$cutpoints[k-1]
          attr(self$cutpoints, "log_Y")[which(self$root$raw_Y == k), 2] = self$cutpoints[k]
        }
        attr(self$cutpoints, "log_Y")[attr(self$root$raw_Y, 'right_cens'), 1] = tail(self$cutpoints, 1)
      } else {
        self$sigma2 = ((sum(self$root$Y^2 * self$root$weights) - 2*sum(self$root$Y*self$root$mu1*self$root$weights) + sum(self$root$mu2 * self$root$weights))) / length(self$root$Y)
      }
      return(invisible(self))
    },

    convergeFit = function(tol = 1e-1, update_sigma2 = FALSE, update_ELBO_progress = TRUE, verbose = TRUE, maxit = Inf) {
      ELBOs = numeric(1000)
      ELBOs[1] = -Inf
      ELBOs[2] = self$root$ELBO
      i = 2
      nodes = Traverse(self$root, 'post-order')
      while (abs(ELBOs[i] - ELBOs[i-1]) > tol & (i-1 <= maxit)) {
        currentInputs = NULL
        for (n in nodes) {
          currentInputs = n$updateFit(currentInputs)
        }

        if (update_sigma2) {
          self$root$updateSigma2()
        }
        i = i+1
        if (all(self$root$sigma2 == 0)) { # if estimate sigma2 is 0, stop
          browser()
          ELBOs[i] = Inf
          break
        }
        if (i > length(ELBOs)) { # double size of ELBOs for efficiency rather than making it bigger each iteration
          ELBOs = c(ELBOs, rep(0, i))
        }
        ELBOs[i] = self$root$ELBO
        if (verbose) {
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

    AddSiblingVEB = function(learner, operator = c("+", "*"), combine_name) { # function to add subtree as sibling to given LEAF node, combining with operator
      # learner is tree to add to self
      # operator is how to combine them
      # name is what to call the combining node

      if (!self$isLeaf) {
        stop("'$AddSiblingVEB' should only be called on a leaf node")
      }

      self_copy = self$clone()

      self$learner = NULL
      self$operator = operator
      self$name = combine_name

      self$AddChildNode(self_copy)
      self$AddChildNode(learner)

      if (all(learner$mu1 == 1 * (operator == "*")) & all(learner$mu2 == 1 * (operator == "*"))) {
        # if adding 'null' node, only need to change the moments stored in our own private field, since everything downstream will be unchanged
        # self_copy$parent$updateMoments()
        # self_copy$updateMoments()
        self_copy$parent$updateMoments()
      } else {
        # otherwise, if adding a real sibling, have to update all ancestors
        # self_copy$parent$updateMomentsAll()
        self_copy$updateMomentsAll()
      }

      return(invisible(self$root)) # return root, so we can assign entire tree to new value (in case we're splitting at root)

    },

    lockSelf = function() { # function node calls on itself to check if it's locked, and lock if necessary
      # if changeToConstant, change fitting function to constant
      if (self$isConstant) {
        self$isLocked = TRUE
        if (self$learner$changeToConstant) {
          self$learner = constLearner
          self$updateFit()
          try({self$updateMomentsAll()}, silent = T) # needed to avoid "attempt to apply non-function" error
        }
      }
      return(invisible(self))
    },

    lockLearners = function() { # lock learners that should be locked
      # self$root$Do(function(node) node$lockSelf(), filterFun = function(node) node$isLeaf & !node$isLocked, traversal ='post-order')
      self$root$Do(function(node) node$lockSelf(), filterFun = function(node) node$isLeaf, traversal ='post-order')

      base_learners = Traverse(self$root, filterFun = function(x) x$isLeaf & !x$isLocked)
      for (base_learner in base_learners) {
        if (base_learner$isRoot || base_learner$parent$isRoot) { # if root or parent is root, not locked, do this to avoid errors in next if statement
          next
        }
        # if only adding or multiplying, and we already added or multiplied with another node, and that node is locked, then lock learner
        if ((base_learner$learner$growMode %in% c("+", "*")) && (base_learner$parent$operator == base_learner$learner$growMode) && base_learner$siblings[[1]]$isLocked) {
          base_learner$isLocked = TRUE
        }
        # if "+*", and already if a "+*" part where both "+" anr "*" parts or locked, then lock learner
        if ((base_learner$learner$growMode == "+*") && (base_learner$parent$operator == "*") && base_learner$siblings[[1]]$isLocked && (base_learner$parent$parent$operator == "+") && base_learner$parent$siblings[[1]]$isLocked) {
          base_learner$isLocked = TRUE
        }
      }

      return(invisible(self$root))
    },

    addLearnerAll = function() { # to each leaf, add a "+" and "*"
      self$root$lockLearners() # lock learners

      base_learners = Traverse(self$root, filterFun = function(x) x$isLeaf & !x$isLocked)
      for (base_learner in base_learners) {
        if (!is.null(base_learner$learner$growMode) && (base_learner$learner$growMode != "NA")) {
          if (base_learner$learner$growMode %in% c("+", "*")) {
            learner_name = paste("mu_", base_learner$root$leafCount, sep = '')
            combine_name = paste("combine_", base_learner$root$leafCount, sep = '')

            add_node = VEBBoostNode$new(learner_name, learner = base_learner$learner)
            add_node$learner$currentFit$mu1 = (add_node$learner$currentFit$mu1 * 0) + (base_learner$learner$growMode == "*" * 1)
            add_node$learner$currentFit$mu2 = (add_node$learner$currentFit$mu2 * 0) + (base_learner$learner$growMode == "*" * 1)
            # add_node$learner$currentFit$KL_div = 0
            add_node$learner$currentFit = list(mu1 = add_node$learner$currentFit$mu1, mu2 = add_node$learner$currentFit$mu2, KL_div = 0)
            base_learner$AddSiblingVEB(add_node, base_learner$learner$growMode, combine_name)
          } else if (base_learner$learner$growMode == "+*") {
            learner_name = paste("mu_", base_learner$root$leafCount, sep = '')
            combine_name = paste("combine_", base_learner$root$leafCount, sep = '')

            add_node = VEBBoostNode$new(learner_name, learner = base_learner$learner)
            add_node$learner$currentFit$mu1 = add_node$learner$currentFit$mu1 * 0
            add_node$learner$currentFit$mu2 = add_node$learner$currentFit$mu2 * 0
            # add_node$learner$currentFit$KL_div = 0
            add_node$learner$currentFit = list(mu1 = add_node$learner$currentFit$mu1, mu2 = add_node$learner$currentFit$mu2, KL_div = 0)
            base_learner$AddSiblingVEB(add_node, "+", combine_name)

            learner_name = paste("mu_", base_learner$root$leafCount, sep = '')
            combine_name = paste("combine_", base_learner$root$leafCount, sep = '')

            mult_node = VEBBoostNode$new(learner_name, learner = add_node$learner)
            mult_node$learner$currentFit$mu1 = mult_node$learner$currentFit$mu1 + 1
            mult_node$learner$currentFit$mu2 = mult_node$learner$currentFit$mu2 + 1
            # mult_node$learner$currentFit$KL_div = 0
            mult_node$learner$currentFit = list(mu1 = mult_node$learner$currentFit$mu1, mu2 = mult_node$learner$currentFit$mu2, KL_div = 0)
            base_learner$children[[1]]$AddSiblingVEB(mult_node, "*", combine_name)
          }
        }
      }

      return(invisible(self$root))
    },

    convergeFitAll = function(tol = 1e-1, update_sigma2 = FALSE, update_ELBO_progress = TRUE, verbose = FALSE, maxit = Inf) {
      self$addLearnerAll()
      self$convergeFit(tol, update_sigma2, update_ELBO_progress, verbose, maxit = maxit)
      return(invisible(self$root))
    },

    predict = function(moment = c(1, 2)) { # function to get prediction on new data
      self$root$Do(function(node) {
        if (1 %in% moment) {
          try({node$pred_mu1 = node$learner$predFunction(node$X_test, node$learner$currentFit, 1)}, silent = TRUE)
        }
        if (2 %in% moment) {
          try({node$pred_mu2 = node$learner$predFunction(node$X_test, node$learner$currentFit, 2)}, silent = TRUE)
        }
      }, filterFun = function(x) x$isLeaf)
      return(invisible(self$root))
    }

  ),
  private = list(
    .Y = NA, # response, only not NA at root
    .sigma2 = NA, # variance, only not NA at root
    .sigma2_prev = NA, # previous variance, only not NA at root
    .mu1 = NA, # current first moment
    .mu2 = NA, # current second moment
    .ELBO_progress = list(-Inf), # list of ELBOs for each iteration of growing and fitting the tree
    .pred_mu1 = NULL, # prediction based on predFunction and given new data (first moment)
    .pred_mu2 = NULL, # prediction based on predFunction and given new data (second moment)
    .isLocked = FALSE, # locked <=> V < V_tol, or both learners directly connected (sibling and parent's sibling) are constant
    .alpha = 0, # used in multi-class learner
    .exposure = 1, # used in poisson.log1pexp
    .ensemble = NULL # used in multi-class learner, either reference to self, or multi-class learner object
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
        mu1 = self$learner$currentFit$mu1
      } else {
        mu1 = private$.mu1
      }
      if (length(mu1) == 1) {
        mu1 = rep(mu1, nrow(as.matrix((self$root$raw_Y))))
      }
      return(mu1)
    },

    mu2 = function(value) {
      if (!missing(value)) {
        stop("`$mu2` cannot be modified directly", call. = FALSE)
      }
      if (self$isLeaf) {
        mu2 = self$learner$currentFit$mu2
      } else {
        mu2 = private$.mu2
      }
      if (length(mu2) == 1) {
        mu2 = rep(mu2, nrow(as.matrix((self$root$raw_Y))))
      }
      return(mu2)
    },

    KL_div = function(value) { # KL divergence from q to g of learner defined by sub-tree with this node as the root
      if (!missing(value)) {
        stop("`$KL_div` cannot be modified directly", call. = FALSE)
      }

      if (self$isLeaf) {
        return(self$learner$currentFit$KL_div)
      }
      return(sum(sapply(self$children, function(x) x$KL_div)))
    },

    X = function(value) { # train predictor for node
      if (missing(value)) {
        if (self$isRoot) {
          return(self$learner$X)
        }
        if (!is.null(self$learner$X)) {
          return(self$learner$X)
        } else {
          return(self$parent$X)
        }
      } else {
        self$learner$currentFit$X = value
      }
    },

    X_test = function(value) { # test predictor for node
      if (missing(value)) {
        if (self$isRoot) {
          return(self$learner$X_test)
        }
        if (!is.null(self$learner$X_test)) {
          return(self$learner$X_test)
        } else {
          return(self$parent$X_test)
        }
      } else {
        self$learner$currentFit$X_test = value
      }
    },

    Y = function(value) { # response for sub-tree
      if (missing(value)) {
        if (self$isRoot) {
          if (self$root$family == "gaussian") {
            return(private$.Y)
          } else if (self$root$family == "binomial") {
            d = self$d
            return((private$.Y - .5 + (self$alpha * d)) / d) # alpha*d needed for multi-class case
          } else if (self$root$family == "negative.binomial") {
            return((private$.Y - private$.exposure) / (2 * self$d * (private$.Y + private$.exposure)))
          } else if (self$root$family == "poisson.log1pexp") {
            return(self$mu1 + (private$.Y - private$.exposure*log1pexp(self$mu1))*(self$root$sigma2) / ((1 + exp(-self$mu1)) * log1pexp(self$mu1)))
          } else if (self$root$family == "aft.loglogistic") {
            d = self$root$d
            s_2d = sqrt(self$root$raw_sigma2)/(2*d)
            res = attr(private$.Y, 'log_Y')[, 1]
            res[attr(private$.Y, 'left_cens')] = attr(private$.Y, 'log_Y')[attr(private$.Y, 'left_cens'), 2] - s_2d[attr(private$.Y, 'left_cens'), 2]
            res[attr(private$.Y, 'right_cens')] = attr(private$.Y, 'log_Y')[attr(private$.Y, 'right_cens'), 1] + s_2d[attr(private$.Y, 'right_cens'), 1]
            res[attr(private$.Y, 'int_cens')] = rowSums((attr(private$.Y, 'log_Y')*d)[attr(private$.Y, 'int_cens'), ]) / rowSums(d[attr(private$.Y, 'int_cens'), ])
            return(res)
          } else if (self$root$family == "ordinal.logistic") {
            d = self$root$d
            rsd = rowSums(d, na.rm = T)
            res = rowSums(d * attr(self$root$cutpoints, "log_Y")) / rsd
            res[which(is.na(attr(self$root$cutpoints, "log_Y")[, 1]))] = self$root$cutpoints[1] - .5/rsd[which(is.na(attr(self$root$cutpoints, "log_Y")[, 1]))]
            res[which(is.na(attr(self$root$cutpoints, "log_Y")[, 2]))] = tail(self$root$cutpoints, 1) + .5/rsd[which(is.na(attr(self$root$cutpoints, "log_Y")[, 2]))]
            return(res)
          } else if (self$root$family == "multinomial.titsias") {
            d = self$ensemble$d
            mu1 = self$ensemble$mu1
            res = mu1[cbind(1:nrow(mu1), attr(self$ensemble$Y, 'which'))] - (.5 / d[, self$root$my_class_index])
            is_my_class = is.na(res) # which obs are in my class
            res[is_my_class] = rowSums(.5 + d[is_my_class, , drop = FALSE]*mu1[is_my_class, , drop = FALSE], na.rm = TRUE) / rowSums(d[is_my_class, , drop = FALSE], na.rm = TRUE)
            return(res)
          } else {
            stop("family must be one of 'gaussian', 'binomial', 'negative.binomial', 'poisson.log1pexp', 'aft.loglogistic', 'ordinal.logistic', or 'multinomial.titsias'")
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
          if (self$family == "gaussian") {
            return(private$.sigma2 / self$weights)
          } else if (self$family == "binomial") {
            return(1 / (self$d * self$weights))
          } else if (self$family == "negative.binomial") {
            return(1 / (self$d * self$weights * (private$.Y + private$.exposure)))
          } else if (self$family == "poisson.log1pexp") {
            return(1 / (((.25*private$.exposure) + .17*private$.Y) * self$weights))
          } else if (self$family == "aft.loglogistic") {
            return(private$.sigma2 / (rowSums(self$d, na.rm = TRUE) * self$weights))
          } else if (self$family == "ordinal.logistic") {
            return(1 / (rowSums(self$d, na.rm = TRUE) * self$weights))
          } else if (self$family == "multinomial.titsias") {
            d = self$ensemble$d
            res = 1 / d[, self$root$my_class_index]
            is_my_class = is.na(res) # which obs are in my class
            res[is_my_class] = 1 / rowSums(d[is_my_class, , drop = FALSE], na.rm = TRUE)
            return(res)
          } else {
            stop("family must be one of 'gaussian', 'binomial', 'negative.binomial', poisson.log1pexp', 'aft.loglogistic', 'ordinal.logistic', or 'multinomial.titsias'")
          }
        }
        if (self$parent$operator == "+") {
          s2 = self$parent$sigma2
          if (length(s2) == 1) {
            s2 = rep(s2, nrow(as.matrix((self$root$raw_Y))))
          }
          return(s2)
        }
        if (self$parent$operator == "*") {
          s2 = self$parent$sigma2 / self$siblings[[1]]$mu2
          if (length(s2) == 1) {
            s2 = rep(s2, nrow(as.matrix((self$root$raw_Y))))
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

    raw_sigma2 = function(value) { # raw value of private$.sigma2 from root
      if (!missing(value)) {
        stop("`$raw_sigma2` cannot be modified directly", call. = FALSE)
      }
      if (self$isRoot) {
        return(private$.sigma2)
      } else {
        return(self$parent$raw_sigma2)
      }
    },

    exposure = function(value) { # exposure variable for poisson or NB, or AFT (for right-censorship info, 1 for censored, 0 for not censored)
      if (missing(value)) {
        if (!(grepl("poisson", self$root$family, ignore.case = TRUE) | (self$root$family == "negative.binomial"))) {
          stop("`$exposure` only used for poisson or negative binomial families")
        }
        if (self$isRoot) {
          e = private$.exposure
          if (length(e) == 1) {
            e = rep(e, length(private$.Y))
          }
          return(e)
        } else {
          return(self$root$exposure)
        }
      } else {
        if (self$isRoot) {
          private$.exposure = value
        } else {
          stop("`$exposure` cannot be modified directly except at the root node", call. = FALSE)
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
      if (self$root$family == "gaussian" || self$root$family == "multinomial.titsias") {
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
          sum(self$root$weights * log(ilogit(xi))) + sum(self$root$weights * (xi / 2) * (d*xi - 1)) + sum(self$root$weights * (self$root$raw_Y - .5) * self$root$mu1) -
            .5*sum(self$root$weights * self$root$mu2 * d) - self$KL_div
        )
      } else if (self$root$family == "negative.binomial") {
        d = self$root$d
        xi = self$root$xi
        return(
          sum(self$root$weights * (self$root$raw_Y + self$root$exposure) * log(ilogit(xi))) + .5*sum(self$root$weights * ((self$root$raw_Y - self$root$exposure) * self$root$mu1 - (self$root$raw_Y + self$root$exposure) * xi)) -
            .5*sum(self$root$weights * d * (self$root$raw_Y + self$root$exposure) * (self$root$mu2 - xi^2)) -
            sum(self$root$weights * (lgamma(self$root$raw_Y + self$root$exposure) - lgamma(self$root$exposure) - lfactorial(self$root$raw_Y))) - self$KL_div
        )
      } else if (self$root$family == "poisson.log1pexp") {
        return(
          -.5*.25*sum(self$root$weights * self$exposure * (self$mu2 - self$mu1^2)) - .5*.17*sum(self$root$weights * self$root$raw_Y * (self$mu2 - self$mu1^2)) +
            sum(self$root$weights * self$root$raw_Y * (log(self$exposure) + loglog1pexp(self$mu1))) - sum(self$root$weights * self$exposure * log1pexp(self$mu1)) -
            sum(self$root$weights * lfactorial(self$root$raw_Y)) - self$KL_div
        )
      } else if (self$root$family == "aft.loglogistic") {
        -.5*log(self$root$raw_sigma2)*sum(self$root$weights[attr(self$root$raw_Y, 'not_cens')]) +
          (.5/sqrt(self$root$raw_sigma2))*(sum(self$root$weights[attr(self$root$raw_Y, 'left_cens')] * (attr(self$root$raw_Y, 'log_Y')[attr(self$root$raw_Y, 'left_cens'), 2] - self$root$mu1[attr(self$root$raw_Y, 'left_cens')])) +
                                             sum(self$root$weights[attr(self$root$raw_Y, 'right_cens')] * (self$root$mu1[attr(self$root$raw_Y, 'right_cens')] - attr(self$root$raw_Y, 'log_Y')[attr(self$root$raw_Y, 'right_cens'), 1])) -
                                             sum(self$root$weights[attr(self$root$raw_Y, 'int_cens')] * attr(self$root$raw_Y, 'log_Y')[attr(self$root$raw_Y, 'int_cens'), ])) -
          .5*sum(cbind(self$root$weights, self$root$weights) * self$root$d * ((cbind(self$root$mu2, self$root$mu2) - 2*cbind(self$root$mu1, self$root$mu1)*attr(self$root$raw_Y, 'log_Y') + attr(self$root$raw_Y, 'log_Y')^2)/self$root$raw_sigma2 - self$xi^2), na.rm = TRUE) +
          sum(self$root$weights[attr(self$root$raw_Y, 'int_cens')] * (attr(self$root$raw_Y, 'log_Y')[attr(self$root$raw_Y, 'int_cens'), 2])/sqrt(self$root$raw_sigma2) + sum(log(-expm1((attr(self$root$raw_Y, 'log_Y')[attr(self$root$raw_Y, 'int_cens'), 1] - attr(self$root$raw_Y, 'log_Y')[attr(self$root$raw_Y, 'int_cens'), 2])/sqrt(self$root$raw_sigma2))))) -
          sum(cbind(self$root$weights, self$root$weights) * (.5*abs(self$root$xi) + cbind(log1pexp(-abs(self$root$xi[, 1])), log1pexp(-abs(self$root$xi[, 2])))), na.rm = TRUE) - self$KL_div
      } else if (self$root$family == "ordinal.logistic") {
        .5*(sum(attr(self$root$cutpoints, 'log_Y')[attr(self$root$raw_Y, 'left_cens'), 2] - self$root$mu1[attr(self$root$raw_Y, 'left_cens')]) +
              sum(self$root$mu1[attr(self$root$raw_Y, 'right_cens')] - attr(self$root$cutpoints, 'log_Y')[attr(self$root$raw_Y, 'right_cens'), 1]) -
              sum(attr(self$root$cutpoints, 'log_Y')[attr(self$root$raw_Y, 'int_cens'), ])) -
          .5*sum(self$root$d*((cbind(self$root$mu2, self$root$mu2) - 2*cbind(self$root$mu1, self$root$mu1)*attr(self$root$cutpoints, 'log_Y') + attr(self$root$cutpoints, 'log_Y')^2) - self$xi^2), na.rm = TRUE) +
          sum(attr(self$root$cutpoints, 'log_Y')[attr(self$root$raw_Y, 'int_cens'), 2]) + sum(log(-expm1((attr(self$root$cutpoints, 'log_Y')[attr(self$root$raw_Y, 'int_cens'), 1] - attr(self$root$cutpoints, 'log_Y')[attr(self$root$raw_Y, 'int_cens'), 2])))) -
          sum(.5*abs(self$root$xi) + cbind(log1pexp(-abs(self$root$xi[, 1])), log1pexp(-abs(self$root$xi[, 2]))), na.rm = TRUE) - self$KL_div
      } else {
        stop("family must be one of 'gaussian', 'binomial', 'negative.binomial', poisson.log1pexp', 'aft.loglogistic', 'ordinal.logistic', or 'multinomial.titsias'")
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
      if (self$root$family == "aft.loglogistic") {
        return(sqrt((1 / self$root$raw_sigma2) * (cbind(self$root$mu2, self$root$mu2) - 2*cbind(self$root$mu1, self$root$mu1)*attr(self$root$raw_Y, 'log_Y') + attr(self$root$raw_Y, 'log_Y')^2)))
      } else if (self$root$family == "ordinal.logistic") {
        return(sqrt((cbind(self$root$mu2, self$root$mu2) - 2*cbind(self$root$mu1, self$root$mu1)*attr(self$root$cutpoints, 'log_Y') + attr(self$root$cutpoints, 'log_Y')^2)))
      }
      return(sqrt(self$root$mu2 + self$alpha^2 - 2*self$alpha*self$root$mu1))
    },

    d = function(value) { # d == 1/xi * (g(xi) - .5), n x K matrix
      if (!missing(value)) {
        stop("'$d' cannot be modified directly", call. = FALSE)
      }
      xi = self$xi
      g_xi = ilogit(xi) # matrix of g(xi_i,k), pre-compute once
      d = ((g_xi - .5) / xi) # matrix of (g(xi_i,k) - .5) / xi_i,k, pre-compute once
      d[xi == 0] = .25 # case of 0/0 (e.g. x_i is all 0), use L'Hopital
      return(d)
    },

    isConstant = function(value) { # uses node's constant check function and current fit to see if it is essentially a constant function
      if (!missing(value)) {
        stop("'$isConstant' cannot be modified directly", call. = FALSE)
      }
      return(self$learner$constCheckFunction(self$learner$currentFit))
    },

    isLocked = function(value) { # locked <=> V < V_tol, or both learners directly connected (sibling and parent's sibling) are constant, or just not growing this learner
      if (missing(value)) {
        if (self$isLeaf) {
          return(private$.isLocked || is.null(self$learner$growMode) || (self$learner$growMode == "NA"))
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
        children_pred_mu1 = lapply(self$children, function(x) x$pred_mu1)
        if (self$operator == "+") {
          return(children_pred_mu1[[1]] + children_pred_mu1[[2]])
        } else {
          return(children_pred_mu1[[1]] * children_pred_mu1[[2]])
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
        children_pred_mu1 = lapply(self$children, function(x) x$pred_mu1)
        children_pred_mu2 = lapply(self$children, function(x) x$pred_mu2)
        return(children_pred_mu2[[1]] + children_pred_mu2[[2]] + (2 * children_pred_mu1[[1]] * children_pred_mu1[[2]]))
      } else {
        children_pred_mu2 = lapply(self$children, function(x) x$pred_mu2)
        return(children_pred_mu2[[1]] * children_pred_mu2[[2]])
      }
    }
  ),
  inherit = data.tree::Node,
  lock_objects = FALSE
)
