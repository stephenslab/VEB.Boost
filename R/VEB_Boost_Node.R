#' @import data.tree
#' @import R6

### Node Object ###

VEBBoostNode <- R6Class(
  "VEBBoostNode",
  public = list(
    operator = NULL, # either "+" or "*" for internal nodes, NULL for terminal nodes
    
    currentFit = NULL, # current fit for fitting function
    
    fitFunction = NULL, # function that takes in predictors X, response Y, variances sigma2, and returns the fit
    # the fit must have fields mu1 (first moment), mu2 (second moment), KL_div (KL divergence from q to g),
    
    predFunction = NULL, # function to predict based on current fit
    
    constCheckFunction = NULL, # function to check if fit is constant
    
    family = "gaussian",
    
    cutpoints = NULL, # cutpoints for ordinal regression
    
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
          self$fitFunction = private$.fitFnConstComp
          self$predFunction = private$.predFnConstComp
          self$constCheckFunction = private$.constCheckFnConstComp
        }
       
        self$currentFit = self$fitFunction(self$X, currentInputs$Y, currentInputs$sigma2, self$currentFit)
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
        a1 = -sum(attr(self$root$raw_Y, 'not_cens'))
        a2 = .5*(sum(attr(self$root$raw_Y, 'log_Y')[attr(self$root$raw_Y, 'left_cens'), 2] - self$root$mu1[attr(self$root$raw_Y, 'left_cens')]) +
                   sum(self$root$mu1[attr(self$root$raw_Y, 'right_cens')] - attr(self$root$raw_Y, 'log_Y')[attr(self$root$raw_Y, 'right_cens'), 1]) -
                   sum(attr(self$root$raw_Y, 'log_Y')[attr(self$root$raw_Y, 'int_cens'), ]))
        a3 = -.5*sum(self$root$d*(cbind(self$root$mu2, self$root$mu2) - 2*cbind(self$root$mu1, self$root$mu1)*attr(self$root$raw_Y, 'log_Y') + attr(self$root$raw_Y, 'log_Y')^2), na.rm = TRUE)
        nll.logscale = function(lS) {
          -(a1*lS + a2*exp(-lS) + a3*exp(-2*lS) +
              sum(attr(self$root$raw_Y, 'log_Y')[attr(self$root$raw_Y, 'int_cens'), 2])/exp(lS) +
              sum(log(-expm1((attr(self$root$raw_Y, 'log_Y')[attr(self$root$raw_Y, 'int_cens'), 1] - attr(self$root$raw_Y, 'log_Y')[attr(self$root$raw_Y, 'int_cens'), 2])/exp(lS)))))
        }
        lS = optim(par = .5*log(self$raw_sigma2), fn = nll.logscale, method = 'Brent', lower = -15, upper = 35)$par
        self$sigma2 = exp(2*lS)
      } else if (self$root$family == "ordinal.logistic") {
        d = self$root$d
        n_k = table(self$root$raw_Y)
        sum_d_k_2 = .5*do.call(rbind, lapply(1:max(self$root$raw_Y), function(k) colSums(d[which(self$root$raw_Y == k), ])))
        sum_dT_k = do.call(rbind, lapply(1:max(self$root$raw_Y), function(k) colSums(sweep(d[which(self$root$raw_Y == k), ], 1, self$root$mu1[which(self$root$raw_Y == k)], "*"))))
        fn = function(theta) {
          res = 0
          for (k in 1:length(theta)) {
            res = res = (theta[k]^2)*(sum_d_k_2[k, 2] + sum_d_k_2[k+1, 1]) - theta[k]*(.5*(n_k[k] - n_k[k+1]) + sum_dT_k[k, 2] + sum_dT_k[k+1, 1])
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
        self$sigma2 = ((sum(self$root$Y^2) - 2*sum(self$root$Y*self$root$mu1) + sum(self$root$mu2))) / length(self$root$Y)
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
      self_copy$X = NULL
      
      self$fitFunction = NULL
      self$currentFit = NULL
      self$predFunction = NULL
      self$constCheckFunction = NULL
      self$operator = operator
      self$name = combine_name
      
      self$AddChildNode(self_copy)
      self$AddChildNode(learner)
      
      if (all(learner$mu1 == 1 * (operator == "*")) & all(learner$mu2 == 1 * (operator == "*"))) {
        # if adding 'null' node, only need to change the moments stored in our own private field, since everything downstream will be unchanged
        self_copy$parent$updateMoments()
      } else {
        # otherwise, if adding a real sibling, have to update all ancestors
        self_copy$parent$updateMomentsAll()
      }
      
      return(invisible(self$root)) # return root, so we can assign entire tree to new value (in case we're splitting at root)
      
    },
    
    lockSelf = function(changeToConstant = TRUE) { # function node calls on itself to check if it's locked, and lock if necessary
      # if changeToConstant, change fitting function to constant
      if (self$isConstant) {
        self$isLocked = TRUE
        if (changeToConstant) {
          self$fitFunction = fitFnConstComp
          self$predFunction = predFnConstComp
          self$constCheckFunction = constCheckFnConstComp
          self$updateFit()
          try({self$updateMomentsAll()}, silent = T) # needed to avoid "attempt to apply non-function" error
        }
      }
      return(invisible(self))
    },
    
    lockLearners = function(growMode = c("+", "*", "+*"), changeToConstant = TRUE) { # lock learners that should be locked
      self$root$Do(function(node) node$lockSelf(changeToConstant), filterFun = function(node) node$isLeaf & !node$isLocked, traversal ='post-order')
      
      base_learners = Traverse(self$root, filterFun = function(x) x$isLeaf & !x$isLocked)
      for (learner in base_learners) {
        if (learner$isRoot || learner$parent$isRoot) { # if root or parent is root, not locked, do this to avoid errors in next if statement
          next
        }
        # if only adding or multiplying, and we already added or multiplied with another node, and that node is locked, then lock learner
        if ((growMode %in% c("+", "*")) && (learner$parent$operator == growMode) && learner$siblings[[1]]$isLocked) {
          learner$isLocked = TRUE
        }
        # if "+*", and already if a "+*" part where both "+" anr "*" parts or locked, then lock learner
        if ((growMode == "+*") && (learner$parent$operator == "*") && learner$siblings[[1]]$isLocked && (learner$parent$parent$operator == "+") && learner$parent$siblings[[1]]$isLocked) {
          learner$isLocked = TRUE
        }
      }
      
      return(invisible(self$root))
    },
    
    unlockLearners = function(changeToConstant = TRUE) { # unlock all (non-constant) learners
      if (changeToConstant) {
        self$root$Do(function(node) node$isLocked = FALSE, filterFun = function(node) node$isLeaf && !node$isConstant && !(node$name %in% c("mu_mrAsh", "mu_RE")))
      } else {
        self$root$Do(function(node) node$isLocked = FALSE, filterFun = function(node) node$isLeaf && !(node$name %in% c("mu_mrAsh", "mu_RE")))
      }
      return(invisible(self$root))
    },
    
    addLearnerAll = function(growMode = c("+", "*", "+*"), changeToConstant = TRUE, unlockLearners = FALSE) { # to each leaf, add a "+" and "*"
      self$root$lockLearners(growMode, changeToConstant) # lock learners
      
      base_learners = Traverse(self$root, filterFun = function(x) x$isLeaf & !x$isLocked)
      for (learner in base_learners) {
        fitFn = learner$fitFunction
        predFn = learner$predFunction
        constCheckFn = learner$constCheckFunction

        if (growMode %in% c("+", "*")) {
          learner_name = paste("mu_", learner$root$leafCount, sep = '')
          combine_name = paste("combine_", learner$root$leafCount, sep = '')
          
          add_fit = list(mu1 = rep(1 * (growMode == "*"), length(learner$Y)), mu2 = rep(1 * (growMode == "*"), length(learner$Y)), KL_div = 0, V = 1)
          add_node = VEBBoostNode$new(learner_name, fitFunction = fitFn, predFunction = predFn, constCheckFunction = constCheckFn, currentFit = add_fit)
          learner$AddSiblingVEB(add_node, growMode, combine_name)
        } else {
          learner_name = paste("mu_", learner$root$leafCount, sep = '')
          combine_name = paste("combine_", learner$root$leafCount, sep = '')
          
          add_fit = list(mu1 = rep(0, length(learner$Y)), mu2 = rep(0, length(learner$Y)), KL_div = 0, V = 1)
          add_node = VEBBoostNode$new(learner_name, fitFunction = fitFn, predFunction = predFn, constCheckFunction = constCheckFn, currentFit = add_fit)
          learner$AddSiblingVEB(add_node, "+", combine_name)
          
          learner_name = paste("mu_", learner$root$leafCount, sep = '')
          combine_name = paste("combine_", learner$root$leafCount, sep = '')
          
          mult_fit = list(mu1 = rep(1, length(learner$Y)), mu2 = rep(1, length(learner$Y)), KL_div = 0, V = 1)
          mult_node = VEBBoostNode$new(learner_name, fitFunction = fitFn, predFunction = predFn, constCheckFunction = constCheckFn, currentFit = mult_fit)
          learner$children[[1]]$AddSiblingVEB(mult_node, "*", combine_name)
        }
      }
      
      if (unlockLearners) {
        self$root$unlockLearners(changeToConstant)
      }
      
      return(invisible(self$root))
    },
    
    convergeFitAll = function(tol = 1e-1, update_sigma2 = FALSE, update_ELBO_progress = TRUE, growMode = "+*", changeToConstant = TRUE, verbose = FALSE, maxit = Inf, unlockLearners = FALSE) {
      self$addLearnerAll(growMode, changeToConstant, unlockLearners)
      self$convergeFit(tol, update_sigma2, update_ELBO_progress, verbose, maxit = maxit)
      return(invisible(self$root))
    },
    
    predict = function(X_new, moment = c(1, 2)) { # function to get prediction on new data
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
        mu1 = self$currentFit$mu1
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
        mu2 = self$currentFit$mu2
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
          } else if (self$root$family == "negative.binomial") {
            return((private$.Y - private$.exposure) / (2 * self$d * (private$.Y + private$.exposure)))
          } else if (self$root$family == "poisson.log1pexp") {
            return(self$mu1 + (private$.Y - private$.exposure*log1pexp(self$mu1))*(self$root$sigma2) / ((1 + exp(-self$mu1)) * log1pexp(self$mu1)))
          } else if (self$root$family == "poisson.exp") {
            # return(self$mu1 + (private$.Y - private$.log1pexp(self$mu1)) / ((1 / self$root$sigma2) * (1 + exp(-self$mu1)) * private$.log1pexp(self$mu1)))
            return(self$mu1 + (private$.Y - private$.exposure*exp(self$mu1))*(self$root$sigma2))
            # xi = self$xi
            # return(xi + self$sigma2*(private$.Y - exp(xi)))
            # xi = self$mu1 + sqrt(1/.05)*sqrt(self$mu2 - self$mu1^2)
            # return(xi - 1 + (private$.Y * self$root$sigma2))
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
          } else {
            stop("family must be one of 'gaussian', 'binomial', 'negative.binomial', 'poisson.log1pexp', 'poisson.exp', 'aft.loglogistic', or 'ordinal.logistic'")
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
          } else if (self$root$family == "negative.binomial") {
            return(1 / (self$d * (private$.Y + private$.exposure)))
          } else if (self$root$family == "poisson.log1pexp") {
            return(1 / ((.25*private$.exposure) + .17*private$.Y))
          } else if (self$root$family == "poisson.exp") {
            # vars = self$root$mu2 - self$root$mu1^2
            # get_sd = function(y, mu1, var) {
            #   if (var <= 0) {
            #     return(1 / (.25 + .17*y))
            #   }
            #   sd = sqrt(var)
            #   h = function(x) -(exp(x)*(exp(x)*y - y*log(1 + exp(x)) + (log(1 + exp(x)))^2)) / ((1 + exp(x))^2 * (log(1 + exp(x)))^2)
            #   h_max = optim(par = list(x = mu1), fn = h, method = 'Brent', lower = mu1 - sqrt(1/.05)*sd, upper = mu1 + sqrt(1/.05)*sd)$value
            #   return(-1/h_max)
            # }
            # return(mapply(get_sd, y = private$.Y, mu1 = self$root$mu1, var = vars, SIMPLIFY = T))
            vars = self$root$mu2 - self$root$mu1^2
            # sigma2 = rep(1, length(vars))
            sigma2 = log(1 + self$root$raw_Y/private$.exposure) + 3
            sigma2[vars > 0] = log(private$.exposure) + self$root$mu1[vars > 0] + sqrt(1/.01 - 1)*sqrt(vars[vars > 0]) # use Cantelli lemma to bound right tail prob <= .01 (since curvature gets worse as x increases)
            # sigma2 = log(1 + self$root$raw_Y) + 3
            return(exp(-sigma2))
          } else if (self$root$family == "aft.loglogistic") {
            return(private$.sigma2 / rowSums(self$root$d, na.rm = TRUE))
          } else if (self$root$family == "ordinal.logistic") {
            return(1 / rowSums(self$root$d, na.rm = TRUE))
          } else {
            stop("family must be one of 'gaussian', 'binomial', 'negative.binomial', poisson.log1pexp', 'poisson.exp', 'aft.loglogistic', or 'ordinal.logistic'")
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
          sum(log(ilogit(xi))) + sum((xi / 2)*(d*xi - 1)) + sum((self$root$raw_Y - .5) * self$root$mu1) -
            .5*sum(self$root$mu2 * d) - self$KL_div
        )
      } else if (self$root$family == "negative.binomial") {
        d = self$root$d
        xi = self$root$xi
        return(
          sum((self$root$raw_Y + self$root$exposure) * log(ilogit(xi))) + .5*sum((self$root$raw_Y - self$root$exposure) * self$root$mu1 - (self$root$raw_Y + self$root$exposure) * xi) -
            .5*sum(d * (self$root$raw_Y + self$root$exposure) * (self$root$mu2 - xi^2)) -
            sum(lgamma(self$root$raw_Y + self$root$exposure) - lgamma(self$root$exposure) - lfactorial(self$root$raw_Y)) - self$KL_div
        )
      } else if (self$root$family == "poisson.log1pexp") {
        return(
          -.5*.25*sum(self$exposure*(self$mu2 - self$mu1^2)) - .5*.17*sum(self$root$raw_Y * (self$mu2 - self$mu1^2)) + sum(self$root$raw_Y * (log(self$exposure) + loglog1pexp(self$mu1))) - sum(self$exposure*log1pexp(self$mu1)) -
            sum(lfactorial(self$root$raw_Y)) - self$KL_div
        )
      } else if (self$root$family == "poisson.exp") {
        # xi = self$xi
        # return(
        #   -.5*sum((self$mu2 - 2*self$mu1*xi + xi^2) / self$sigma2) + sum((private$.Y - exp(xi))*(self$mu1 - xi)) + sum(self$root$raw_Y * xi) - sum(exp(xi)) -
        #     self$KL_div
        # )
        # return(
        #   -.5*sum((self$mu2 - self$mu1^2) / self$sigma2) + sum(self$root$raw_Y * private$.loglog1pexp(self$mu1)) - sum(exp(private$.loglog1pexp(self$mu1))) -
        #     self$KL_div
        # )
        return(
          -.5*sum((self$mu2 - self$mu1^2) / self$sigma2) + sum(self$root$raw_Y * (log(self$exposure) + self$mu1)) - sum(self$exposure*exp(self$mu1)) -
            sum(lfactorial(self$root$raw_Y)) - self$KL_div
        )
      } else if (self$root$family == "aft.loglogistic") {
        -.5*log(self$root$raw_sigma2)*sum(attr(self$root$raw_Y, 'not_cens')) +
          (.5/sqrt(self$root$raw_sigma2))*(sum(attr(self$root$raw_Y, 'log_Y')[attr(self$root$raw_Y, 'left_cens'), 2] - self$root$mu1[attr(self$root$raw_Y, 'left_cens')]) +
                                             sum(self$root$mu1[attr(self$root$raw_Y, 'right_cens')] - attr(self$root$raw_Y, 'log_Y')[attr(self$root$raw_Y, 'right_cens'), 1]) -
                                             sum(attr(self$root$raw_Y, 'log_Y')[attr(self$root$raw_Y, 'int_cens'), ])) -
          .5*sum(self$root$d*((cbind(self$root$mu2, self$root$mu2) - 2*cbind(self$root$mu1, self$root$mu1)*attr(self$root$raw_Y, 'log_Y') + attr(self$root$raw_Y, 'log_Y')^2)/self$root$raw_sigma2 - self$xi^2), na.rm = TRUE) +
          sum(attr(self$root$raw_Y, 'log_Y')[attr(self$root$raw_Y, 'int_cens'), 2])/sqrt(self$root$raw_sigma2) + sum(log(-expm1((attr(self$root$raw_Y, 'log_Y')[attr(self$root$raw_Y, 'int_cens'), 1] - attr(self$root$raw_Y, 'log_Y')[attr(self$root$raw_Y, 'int_cens'), 2])/sqrt(self$root$raw_sigma2)))) -
          sum(.5*abs(self$root$xi) + cbind(log1pexp(-abs(self$root$xi[, 1])), log1pexp(-abs(self$root$xi[, 2]))), na.rm = TRUE) - self$KL_div
      } else if (self$root$family == "ordinal.logistic") {
        .5*(sum(attr(self$root$cutpoints, 'log_Y')[attr(self$root$raw_Y, 'left_cens'), 2] - self$root$mu1[attr(self$root$raw_Y, 'left_cens')]) +
              sum(self$root$mu1[attr(self$root$raw_Y, 'right_cens')] - attr(self$root$cutpoints, 'log_Y')[attr(self$root$raw_Y, 'right_cens'), 1]) -
              sum(attr(self$root$cutpoints, 'log_Y')[attr(self$root$raw_Y, 'int_cens'), ])) -
          .5*sum(self$root$d*((cbind(self$root$mu2, self$root$mu2) - 2*cbind(self$root$mu1, self$root$mu1)*attr(self$root$cutpoints, 'log_Y') + attr(self$root$cutpoints, 'log_Y')^2) - self$xi^2), na.rm = TRUE) +
          sum(attr(self$root$cutpoints, 'log_Y')[attr(self$root$raw_Y, 'int_cens'), 2]) + sum(log(-expm1((attr(self$root$cutpoints, 'log_Y')[attr(self$root$raw_Y, 'int_cens'), 1] - attr(self$root$cutpoints, 'log_Y')[attr(self$root$raw_Y, 'int_cens'), 2])))) -
          sum(.5*abs(self$root$xi) + cbind(log1pexp(-abs(self$root$xi[, 1])), log1pexp(-abs(self$root$xi[, 2]))), na.rm = TRUE) - self$KL_div
      } else {
        stop("family must be one of 'gaussian', 'binomial', 'negative.binomial', poisson.log1pexp', 'poisson.exp', 'aft.loglogistic', or 'ordinal.logistic'")
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
    
    # xi = function(value) { # optimal variational parameters, set to +sqrt(mu2)
    #   if (!missing(value)) {
    #     stop("`$xi` cannot be modified directly", call. = FALSE)
    #   }
    #   if (self$root$family == 'binomial') {
    #     return(sqrt(self$root$mu2 + self$alpha^2 - 2*self$alpha*self$root$mu1))
    #   } else if (self$root$family == 'poisson') {
    #     f = function(xi, mu=1, sigma2=1, k=sqrt(20)) {
    #       -(.5*exp(mu + k*sqrt(sigma2))*(sigma2 + mu^2) + (exp(xi) - xi*exp(mu + k*sqrt(sigma2)))*mu + (exp(xi) - xi*exp(xi) + .5*exp(mu + k*sqrt(sigma2))*xi^2))
    #     }
    #     get_xi = function(mu, sigma2) {
    #       if (sigma2 <= 0) {
    #         return(mu)
    #       }
    #       optim(list(xi = mu), fn = function(x) f(x, mu, sigma2), method = 'Brent', lower = mu - sqrt(20 * sigma2), upper = mu + sqrt(20 * sigma2))$par
    #     }
    #     mu = self$mu1
    #     sigma2 = self$mu2 - self$mu1^2
    #     xis = mapply(get_xi, mu = mu, sigma2 = sigma2, SIMPLIFY = T)
    #     return(xis)
    #   }
    # },
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
