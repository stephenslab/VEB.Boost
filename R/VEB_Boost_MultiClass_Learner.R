#' @import data.tree
#' @importFrom parallel mclapply

### Multi-Class Learner Object ###

### This class deals with multi-class learners. There's a tree associated with each of the K classes. A quadratic bound is used on the
### soft-max function, which results in a lower bound for the ELBO

VEBBoostMultiClassLearner <- R6::R6Class(
  "VEBBoostMultiClassLearner",
  public = list(
    classLearners = list(), # list of K learners, each is a VEBBoostNode object

    family = "multinomial.bouchard", # either "multinomial.bouchard" or "multinomial.titsias"

    classes = NULL,

    alpha = 0,

    updateAlpha = function() {
      # make n x K matrix of d_i,k
      d = self$d # store

      K = ncol(d)
      m = sapply(self$classLearners, function(x) x$root$mu1) # n x K matrix of first moments

      self$alpha = ((K/2 - 1) + rowSums(d * m)) / rowSums(d)
      return(invisible(self))
    },

    addclassLearners = function(learnerList) { # this function adds the learners from learnerList to the list of learners (and sets the ensemble to self)
      # set `$ensemble` to self for each learner
      for (learner in learnerList) {
        learner$ensemble = self
      }
      self$classLearners = c(self$classLearners, learnerList)
      return(invisible(self))
    },

    convergeFit = function(tol = 1e-1, update_ELBO_progress = TRUE, verbose = FALSE) {
      ELBOs = numeric(1000)
      ELBOs[1] = -Inf
      ELBOs[2] = self$ELBO
      i = 2
      while (abs(ELBOs[i] - ELBOs[i-1]) > tol) {
        for (classLarner in self$classLearners) { # update each node of each learner, could be parallelized and update alpha after updating all learners
          classLarner$root$convergeFit(tol = tol, update_sigma2 = FALSE, update_ELBO_progress = FALSE, verbose = FALSE, maxit = 1)
          # self$updateAlpha() # update alpha after updating learner (could do this just once at the end, too)
        }
        # self$classLearners = foreach(learner = self$classLearners, .packages = c('bigmemory', 'data.tree'), .combine = 'c') %dopar% { # update each node of each learner, could be parallelized and update alpha after updating all learners
        #   learner$root$Do(function(x) try({x$updateFit()}, silent = T), traversal = 'post-order')
        #   learner
        # }
        # mclapply(self$classLearners, private$.updateLearnerFn, mc.cores = self$mc.cores)
        # have to re-point ensemble to self, since the environments get messed up after the parallel foreach call
        # for (learner in self$classLearners) {
        #   learner$ensemble = self
        # }
        if (self$family == "multinomial.bouchard") {
          self$updateAlpha()
        }

        i = i + 1
        if (i > length(ELBOs)) { # double size of ELBOs for efficiency rather than making it bigger each iteration
          ELBOs = c(ELBOs, rep(0, i))
        }
        ELBOs[i] = self$ELBO
        if (verbose) {
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

    convergeFitAll = function(tol = 1e-1, update_ELBO_progress = TRUE, verbose = FALSE) {
      for (classLearner in self$classLearners) {
        # classLearner = classLearner$addLearnerAll()
        classLearner$addLearnerAll()
      }
      self$convergeFit(tol, update_ELBO_progress, verbose)
      return(invisible(self))
    },

    predict = function(moment = c(1, 2)) { # function to get prediction on new data
      for (classLearner in self$classLearners) {
        classLearner$root$predict(moment)
      }
      return(invisible(self))
    }

  ),
  private = list(
    .Y = NA, # response
    .mu1 = NA, # current first moment
    .mu2 = NA, # current second moment
    .ELBO_progress = list(-Inf), # list of ELBOs for each iteration of growing and fitting the tree
    .pred_mu1 = NULL, # prediction based on predFunction and given new data (first moment)
    .pred_mu2 = NULL # prediction based on predFunction and given new data (second moment)
  ),

  active = list(

    Y = function(value) {
      if (missing(value)) {
        return(private$.Y)
      } else {
        private$.Y = value
      }
    },

    mu1 = function(value) {
      if (!missing(value)) {
        stop("`$mu1` cannot be modified directly", call. = FALSE)
      }
      return(sapply(self$classLearners, function(x) x$root$mu1))
    },

    mu2 = function(value) {
      if (!missing(value)) {
        stop("`$mu2` cannot be modified directly", call. = FALSE)
      }
      return(sapply(self$classLearners, function(x) x$root$mu2))
    },

    KL_div = function(value) { # KL divergence from q to g of learners
      if (!missing(value)) {
        stop("`$KL_div` cannot be modified directly", call. = FALSE)
      }

      return(sum(sapply(self$classLearners, function(x) x$root$KL_div)))
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
    #   K = length(self$classLearners)
    #   d = self$d
    #   a = self$alpha
    #   if (length(a) == 1) {
    #     a = rep(a, length(self$Y))
    #   }
    #   b = .5 - (a * d)
    #   ELBO = 0
    #   for (learner in self$classLearners) {
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
      if (self$family == "multinomial.bouchard") {
        K = length(self$classLearners)
        d = self$d
        a = self$alpha
        if (length(a) == 1) {
          a = rep(a,  length(self$Y))
        }
        b = .5 - (a * d)
        ELBO = 0
        for (classLearner in self$classLearners) {
          ELBO = ELBO + sum(classLearner$root$mu1 * classLearner$root$raw_Y)
        }
        xi = self$xi
        c = ((1 - K/2) * a) - (rowSums(xi) / 2) + (rowSums(d * (a^2 - xi^2)) / 2) + rowSums(log(1 + exp(xi)))
        # c = ((1 - K/2) * a) - (rowSums(xi) / 2) + (rowSums(d * (a^2 - xi^2)) / 2) + rowSums(apply(xi, 2, log1pexp))
        m1 = self$mu1
        m2 = self$mu2
        ELBO = ELBO - (sum(m2 * d) / 2) - sum(m1 * b) - sum(c)
        ELBO = ELBO - self$KL_div
        return(ELBO)
      } else if (self$family == "multinomial.titsias") {
        K = length(self$classLearners)
        xi = self$xi
        mu1 = self$mu1
        mu2 = self$mu2
        ELBO = .5 * (K*sum(mu1[cbind(1:nrow(mu1), attr(self$Y, 'which'))]) - sum(mu1)) - sum(apply(-xi, 2, log1pexp), na.rm = TRUE) - .5*sum(xi, na.rm = TRUE)
        ELBO = ELBO - self$KL_div
        return(ELBO)
      } else {
        stop("`$family` must be either 'multinomial.bouchard' or 'multinomial.titsias")
      }
    },

    xi = function(value) { # optimal variational parameters, set to +sqrt(mu2)
      if (!missing(value)) {
        stop("`$xi` cannot be modified directly", call. = FALSE)
      }
      #xi = sapply(self$classLearners, function(x) sqrt(x$root$mu2 + x$root$alpha^2 - 2*x$root$alpha*x$root$mu1))
      if (self$family == "multinomial.bouchard") {
        return(sapply(self$classLearners, function(x) x$xi))
      } else if (self$family == "multinomial.titsias") {
        mu1 = self$mu1
        mu2 = self$mu2
        # xi = do.call(rbind, lapply(1:nrow(mu1), function(i) {
        #   res = mu2[i, attr(self$Y, 'which')[i]] - 2*mu1[i, attr(self$Y, 'which')[i]]*mu1[i, ] + mu2[i, ]
        #   res[attr(self$Y, 'which')[i]] = NA
        #   return(res)
        # }))
        xi = matrix(mu2[cbind(1:nrow(mu2), attr(self$Y, 'which'))], nrow = nrow(mu2), ncol = ncol(mu2), byrow = FALSE) -
          2*matrix(mu1[cbind(1:nrow(mu2), attr(self$Y, 'which'))], nrow = nrow(mu1), ncol = ncol(mu1), byrow = FALSE)*mu1 + mu2
        xi[cbind(1:nrow(xi), attr(self$Y, 'which'))] = NA
        return(sqrt(xi))
      } else {
        stop("`$family` must be either 'multinomial.bouchard' or 'multinomial.titsias")
      }
    },

    d = function(value) { # d == 1/xi * (g(xi) - .5), n x K matrix
      if (!missing(value)) {
        stop("'$d' cannot be modified directly", call. = FALSE)
      }
      # g_xi = g(self$xi) # matrix of g(xi_i,k), pre-compute once
      # d = ((g_xi - .5) / self$xi) # matrix of (g(xi_i,k) - .5) / xi_i,k, pre-compute once
      # d[self$xi == 0] = .25 # case of 0/0 (e.g. x_i is all 0), use L'Hopital
      if (self$family == "multinomial.bouchard") {
        return(sapply(self$classLearners, function(x) x$d))
      } else if (self$family == "multinomial.titsias") {
        xi = self$xi
        g_xi = ilogit(xi) # matrix of g(xi_i,k), pre-compute once
        d = ((g_xi - .5) / xi) # matrix of (g(xi_i,k) - .5) / xi_i,k, pre-compute once
        d[xi == 0] = .25 # case of 0/0 (e.g. x_i is all 0), use L'Hopital
        return(d)
      } else {
        stop("`$family` must be either 'multinomial.bouchard' or 'multinomial.titsias")
      }
    },

    pred_mu1 = function(value) { # predicted first moment given new data
      if (!missing(value)) {
        stop("`$pred_mu1` cannot be modified directly except at leaf nodes by calling `$predict.veb`", call. = FALSE)
      }
      sapply(self$classLearners, function(x) x$root$pred_mu1)
    },

    pred_mu2 = function(value) { # predicted second moment given new data
      if (!missing(value)) {
        stop("`$pred_mu2` cannot be modified directly except at leaf nodes by calling `$predict.veb`", call. = FALSE)
      }
      sapply(self$classLearners, function(x) x$root$pred_mu2)
    }
  ),
  lock_objects = FALSE
)
