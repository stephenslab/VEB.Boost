#' @import data.tree
#' @importFrom parallel mclapply

### Multi-Class Learner Object ###

### This class deals with multi-class learners. There's a tree associated with each of the K classes. A quadratic bound is used on the
### soft-max function, which results in a lower bound for the ELBO

VEBBoostMultiClassLearner <- R6::R6Class(
  "VEBBoostMultiClassLearner",
  public = list(
    learners = list(), # list of K learners, each is a VEBBoostNode object

    X = NULL,

    Y = NULL,

    classes = NULL,

    mc.cores = 1, # parallel cores to use

    alpha = 0,

    updateAlpha = function() {
      # make n x K matrix of d_i,k
      d = self$d # store

      K = ncol(d)
      m = sapply(self$learners, function(x) x$root$mu1) # n x K matrix of first moments

      self$alpha = ((K/2 - 1) + rowSums(d * m)) / rowSums(d)
      return(invisible(self))
    },

    addLearners = function(learnerList) { # this function adds the learners from learnerList to the list of learners (and sets the ensemble to self)
      # set `$ensemble` to self for each learner
      for (learner in learnerList) {
        learner$ensemble = self
      }
      self$learners = c(self$learners, learnerList)
      return(invisible(self))
    },

    convergeFit = function(tol = 1e-3, update_ELBO_progress = TRUE, verbose = FALSE) {
      ELBOs = numeric(1000)
      ELBOs[1] = -Inf
      ELBOs[2] = self$ELBO
      i = 2
      while (abs(ELBOs[i] - ELBOs[i-1]) > tol) {
        # for (learner in self$learners) { # update each node of each learner, could be parallelized and update alpha after updating all learners
        #   learner$root$Do(function(x) try({x$updateFit()}, silent = T), traversal = 'post-order')
        #   self$updateAlpha() # update alpha after updating learner
        # }
        # self$learners = foreach(learner = self$learners, .packages = c('bigmemory', 'data.tree'), .combine = 'c') %dopar% { # update each node of each learner, could be parallelized and update alpha after updating all learners
        #   learner$root$Do(function(x) try({x$updateFit()}, silent = T), traversal = 'post-order')
        #   learner
        # }
        mclapply(self$learners, private$.updateLearnerFn, mc.cores = self$mc.cores)
        # have to re-point ensemble to self, since the environments get messed up after the parallel foreach call
        # for (learner in self$learners) {
        #   learner$ensemble = self
        # }
        self$updateAlpha()

        i = i+1
        if (i > length(ELBOs)) { # double size of ELBOs for efficiency rather than making it bigger each iteration
          ELBOs = c(ELBOs, rep(0, i))
        }
        ELBOs[i] = self$ELBO
        if (verbose & ((i %% 1) == 0)) {
          cat(paste("iteration: ", i-2, "\t ELBO: ", ELBOs[i], "\t ELBO change: ", ELBOs[i] - ELBOs[i-1], sep = ""))
          cat("\n")
        }
      }
      ELBOs = ELBOs[2:i]

      if (update_ELBO_progress) {
        self$ELBO_progress = ELBOs
      }

      return(invisible(self))
    },

    convergeFitAll = function(tol = 1e-3, update_ELBO_progress = TRUE, verbose = FALSE) {
      # for (learner in self$learners) {
      #   learner = learner$addLearnerAll()
      # }
      mclapply(self$learners, private$.addLearnerFn, changeToConstant = changeToConstant, mc.cores = self$mc.cores)
      self$convergeFit(tol, update_ELBO_progress, verbose)
      return(invisible(self))
    },

    predict.veb = function(X_new, moment = c(1, 2)) { # function to get prediction on new data
      for (learner in self$learners) {
        learner$root$Do(function(node) {
          if (1 %in% moment) {
            node$pred_mu1 = node$predFunction(X_new, node$currentFit, 1)
          }
          if (2 %in% moment) {
            node$pred_mu2 = node$predFunction(X_new, node$currentFit, 2)
          }
        }, filterFun = function(x) x$isLeaf)
      }
      return(invisible(self))
    }

  ),
  private = list(
    .mu1 = NA, # current first moment
    .mu2 = NA, # current second moment
    .ELBO_progress = list(-Inf), # list of ELBOs for each iteration of growing and fitting the tree
    .pred_mu1 = NULL, # prediction based on predFunction and given new data (first moment)
    .pred_mu2 = NULL, # prediction based on predFunction and given new data (second moment)
    .updateLearnerFn = function(learner) { # function for learner to update fit
      learner$root$Do(function(x) x$updateFit(), traversal = 'post-order')
      return(learner)
    },
    .addLearnerFn = function(learner, changeToConstant) { # function for learner to add nodes
      learner$addLearnerAll(changeToConstant)
      return(learner)
    }
  ),

  active = list(

    mu1 = function(value) {
      if (!missing(value)) {
        stop("`$mu1` cannot be modified directly", call. = FALSE)
      }
      return(sapply(self$learners, function(x) x$root$mu1))
    },

    mu2 = function(value) {
      if (!missing(value)) {
        stop("`$mu2` cannot be modified directly", call. = FALSE)
      }
      return(sapply(self$learners, function(x) x$root$mu2))
    },

    KL_div = function(value) { # KL divergence from q to g of learners
      if (!missing(value)) {
        stop("`$KL_div` cannot be modified directly", call. = FALSE)
      }

      return(sum(sapply(self$learners, function(x) x$root$KL_div)))
    },

    ELBO_progress = function(value) { # ELBO progress of tree
      if (missing(value)) {
        return(private$.ELBO_progress)
      } else {
        private$.ELBO_progress[[length(private$.ELBO_progress) + 1]] = value
      }
    },

    # ELBO = function(value) { # ELBO for entire tree
    #   if (!missing(value)) {
    #     stop("`$ELBO` cannot be modified directly", call. = FALSE)
    #   }
    #   K = length(self$learners)
    #   d = self$d
    #   a = self$alpha
    #   if (length(a) == 1) {
    #     a = rep(a, length(self$Y))
    #   }
    #   b = .5 - (a * d)
    #   ELBO = 0
    #   for (learner in self$learners) {
    #     ELBO = ELBO + sum(learner$root$mu1 * learner$root$raw_Y)
    #   }
    #   xi = self$xi
    #   m1 = self$mu1
    #   m2 = self$mu2
    #   for (i in 1:length(self$Y)) {
    #     d_i = d[i, ]
    #     b_i = b[i, ]
    #     xi_i = xi[i, ]
    #     a_i = a[i]
    #     c_i = ((1 - K/2) * a_i) - (sum(xi_i) / 2) + (sum(d_i * (a_i^2 - xi_i^2)) / 2) + sum(log(1 + exp(xi_i)))
    #     m1_i = m1[i, ] # first moments for K learners for this ob
    #     m2_i = m2[i, ] # second moments for K learners for this ob
    #     ELBO = ELBO - (sum(m2_i * d_i) / 2) - sum(m1_i * b_i) - c_i
    #   }
    #   ELBO = ELBO - self$KL_div
    #   return(ELBO)
    # },

    ELBO = function(value) { # ELBO for entire tree
      if (!missing(value)) {
        stop("`$ELBO` cannot be modified directly", call. = FALSE)
      }
      K = length(self$learners)
      d = self$d
      a = self$alpha
      if (length(a) == 1) {
        a = rep(a,  length(self$Y))
      }
      b = .5 - (a * d)
      ELBO = 0
      for (learner in self$learners) {
        ELBO = ELBO + sum(learner$root$mu1 * learner$root$raw_Y)
      }
      xi = self$xi
      c = ((1 - K/2) * a) - (rowSums(xi) / 2) + (rowSums(d * (a^2 - xi^2)) / 2) + rowSums(log(1 + exp(xi)))
      m1 = self$mu1
      m2 = self$mu2
      ELBO = ELBO - (sum(m2 * d) / 2) - sum(m1 * b) - sum(c)
      ELBO = ELBO - self$KL_div
      return(ELBO)
    },

    xi = function(value) { # optimal variational parameters, set to +sqrt(mu2)
      if (!missing(value)) {
        stop("`$xi` cannot be modified directly", call. = FALSE)
      }
      #xi = sapply(self$learners, function(x) sqrt(x$root$mu2 + x$root$alpha^2 - 2*x$root$alpha*x$root$mu1))
      return(sapply(self$learners, function(x) x$xi))
    },

    d = function(value) { # d == 1/xi * (g(xi) - .5), n x K matrix
      if (!missing(value)) {
        stop("'$d' cannot be modified directly", call. = FALSE)
      }
      # g_xi = g(self$xi) # matrix of g(xi_i,k), pre-compute once
      # d = ((g_xi - .5) / self$xi) # matrix of (g(xi_i,k) - .5) / xi_i,k, pre-compute once
      # d[self$xi == 0] = .25 # case of 0/0 (e.g. x_i is all 0), use L'Hopital
      return(sapply(self$learners, function(x) x$d))
    },

    pred_mu1 = function(value) { # predicted first moment given new data
      if (!missing(value)) {
        stop("`$pred_mu1` cannot be modified directly except at leaf nodes by calling `$predict.veb`", call. = FALSE)
      }
      sapply(self$learners, function(x) x$root$pred_mu1)
    },

    pred_mu2 = function(value) { # predicted second moment given new data
      if (!missing(value)) {
        stop("`$pred_mu2` cannot be modified directly except at leaf nodes by calling `$predict.veb`", call. = FALSE)
      }
      sapply(self$learners, function(x) x$root$pred_mu2)
    }
  ),
  lock_objects = FALSE
)
