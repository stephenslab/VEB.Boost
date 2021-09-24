#' set up a general trend filtering matrix
#' @param t vector of length n specifying locations of data points on x axis
#' @param br vector of length (p-1) specifying break points on x axis (ie where changepoints can occur)
#' By default br=t which allows breaks to occur between every data point. Note that internally duplicate elements of br are removed.
#' @param order non-negative integer indicating order of trend filtering basis (0 is changepoint basis and is the only case we test and use)
#' @importFrom Matrix sparseVector
#' @importFrom matrixStats binCounts
#' @keywords internal
make_tfg_matrix = function(t, br = t, order = 0) {
  br = unique(sort(br))
  n = length(t) # number of data points
  p = length(br) + 1 # number of bins specified by breaks
  X <- numeric(0)
  attr(X, "nrow") <- n
  attr(X, "ncol") <- p
  attr(X, "matrix.type") = "tfg_matrix"
  attr(X, "order") = order
  attr(X,"br") <- br
  order_t <- order(t)
  counts = binCounts(t, bx = c(-.Machine$double.xmax, br, .Machine$double.xmax), right=TRUE) #hist(t, breaks = c(-Inf,br,Inf), plot=FALSE)$counts
  attr(X, "null_bin") <- which.max(counts)
  attr(X,"bin_to_t") <- cumsum(counts)
  t_to_bin <- .bincode(t,breaks = c(-Inf,br,Inf))
  i = which(t_to_bin != attr(X, "null_bin"))
  nz = t_to_bin[i]
  attr(X,"t_to_bin") <- sparseVector(nz, i, length = n) # .bincode for t as a sparse vector, 0 signifies the "null bin"
  if (length(attr(X,"t_to_bin")@i) / attr(X,"t_to_bin")@length > .5) { # if more efficient to store as regular vector, change to regular vector
    attr(X,"t_to_bin") = rep(attr(X, "null_bin"), n)
    attr(X,"t_to_bin")[i] = nz
  }
  if ((attr(X, "null_bin") == 1) | (attr(X, "null_bin") == 2 & diff(head(attr(X, "bin_to_t"), 2)) == 0)) { # if biggest group is first bin
    attr(X, "order_t_low") = numeric()
  } else {
    attr(X, "order_t_low") = order_t[1:(attr(X,"bin_to_t")[attr(X, "null_bin")-1])] # ordering of obs less than biggest group
  }
  if ((attr(X, "null_bin") == p) | (attr(X, "null_bin") == p-1 & diff(tail(attr(X, "bin_to_t"), 2)) == 0)) { # if biggest group is last bin
    attr(X, "order_t_high") = numeric()
  } else {
    attr(X, "order_t_high") = order_t[(attr(X, "bin_to_t")[attr(X, "null_bin")]+1):n] # ordering of obs more than biggest group
  }
  # attr(X,"bin_to_t") <- cumsum(matrixStats::binCounts(t, bx = c(-Inf,br,Inf), right=TRUE)) # might be faster, but have to deal w/ Inf
  attr(X,"scaled:center") <- 0
  attr(X,"scaled:scale") <- 1
  return(X)
}

is.tfg_matrix=function(X) {
  ifelse(is.null(attr(X, "matrix.type")), FALSE, attr(X, "matrix.type") == "tfg_matrix")
}

# over-write make_stumps_matrix to not rely on susieR::
# also change Xtrain to be a list (allows for different lengths of breaks)
# include_linear now supports a logical vector input with the same length as ncol(X)
# for those entries that are TRUE, it includes those variables as linear terms
make_stumps_matrix = function(X, include_linear, include_stumps, Xtrain = NULL) {
  if (is.null(Xtrain)) {
    Xtrain = lapply(1:ncol(X), function(j) X[, j])
  }
  
  if (length(include_linear) == 1) { # change include_linear to be a logical vector
    include_linear = rep(include_linear, ncol(X))
  }
  
  xl = list() # initialize
  if(any(include_linear)){ # include X as a regular matrix first
    X_linear = X[, include_linear]
    attr(X_linear, "nrow") <- nrow(X_linear)
    attr(X_linear, "ncol") <- ncol(X_linear)
    attr(X_linear, "scaled:center") <- rep(0, ncol(X_linear))
    attr(X_linear, "scaled:scale") <- rep(1, ncol(X_linear))
    X_linear2 = X_linear^2
    attr(X_linear, "X2") <- X_linear2
    xl = c(xl, list(X_linear))
  }
  
  if (length(include_stumps) == 1) { # change include_stumps to be a logical vector
    include_stumps = rep(include_stumps, ncol(X))
  }
  
  if (any(include_stumps)) {
    for(i in 1:ncol(X)){
      if (include_stumps[i]) {
        xl = c(xl, list(make_tfg_matrix(X[, i], Xtrain[[i]])))
      }
    }
  }
  
  return(xl)
}

# over-write stumps multiplication (WHY IS THERE A -1*.... ?!?!?)
#' @title Compute unscaled X \%*\% b using the special structure of trend filtering
#' @param X a tfg_matrix created by make_tfg_matrix
#' @param b a p vector of the changes at each change point
#' @return an n vector of the means at each data point
#' @importFrom spatstat.utils revcumsum
#' @keywords internal
compute_tfg_Xb = function(X,b){
  order = get_order(X)
  for(i in 1:(order+1)){
    #b = rev(cumsum(rev(b))) # computes mean in each bin
    b = revcumsum(b) # faster than rev(cumsum(rev(b)))
  }
  if (grepl('sparseVector', class(attr(X, "t_to_bin")), ignore.case = T)) { # if using sparse representation
    Xb = rep(b[attr(X, "null_bin")], attr(X, "nrow")) # initialize with null bin
    Xb[attr(attr(X, "t_to_bin"), 'i')] = b[attr(attr(X, "t_to_bin"), 'x')] # fill in non-null bin values
  } else {
    Xb = b[attr(X,"t_to_bin")]
  }
  return(Xb) #  maps bin means to a mean for each datapoint
}

#' @title Compute t(X) \%*\% y using the special structure of trend filtering
#' @param X a tfg_matrix created by make_tfg_matrix
#' @param y an n vector of data
#' @return a p vector
#' @keywords internal
compute_tfg_Xty = function(X,y){
  order = get_order(X)
  if (grepl('sparseVector', class(attr(X, "t_to_bin")), ignore.case = T)) { # if using sparse representation
    y = c(y[attr(X, "order_t_low")], y[-attr(attr(X, "t_to_bin"), "i")], y[attr(X, "order_t_high")])
  } else {
    y = y[c(attr(X, "order_t_low"), which(attr(X, "t_to_bin") == attr(X, "null_bin")), attr(X, "order_t_high"))]
  }
  for (i in 1:(order+1)){
    y = cumsum(y)
  }
  y = y[attr(X,"bin_to_t")]
  # corner case where first bucket has no obs
  if (attr(X, "bin_to_t")[1] == 0) {
    y = c(0, y)
  }
  return(y)
}