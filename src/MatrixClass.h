#ifndef INCLUDE_MatrixClass_h
#define INCLUDE_MatrixClass_h

#include <RcppArmadillo.h>
#include <memory>
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// Matrix block abstract class, i.e. an element of X_stumps (either numeric matrix, or tfg matrix)
class MatrixBlock {
  public: 
    MatrixBlock() {}
    
    ~MatrixBlock() {}
    
    // method to compute (X - X_avg) * b
    virtual vec compute_Xb(const vec& b, const vec& X_avg) = 0;
    
    // method to compute (X - X_avg)' * y
    virtual vec compute_Xty(const vec& y, const vec& X_avg) = 0;
    
    // method to compute (X - X_avg)^2 * b
    virtual vec compute_X2b(const vec& b, const vec& X_avg) = 0;
    
    // method to compute (X' - X_avg')^2 * y
    virtual vec compute_X2ty(const vec& y, const vec& X_avg) = 0;
    
    unsigned int nrow, ncol;
    unsigned int col_offset; // column offset in entire matrix, i.e. sum of # columns before this block
};

// Numeric matrix block, i.e. first element of X_stumps when including linear
class NumericMatrixBlock: public MatrixBlock {
  public:
    ~NumericMatrixBlock() {}
    
    vec compute_Xb(const vec& b, const vec& X_avg) {
      vec res = X*b - as_scalar(dot(b, X_avg));
      return(res);
    }
    
    vec compute_Xty(const vec& y, const vec& X_avg) {
      vec res = trans((y.t() * X)) - (X_avg  * sum(y));
      return(res);
    }
    
    vec compute_X2b(const vec& b, const vec& X_avg) {
      vec X2b(X.n_rows); X2b = X2*b;
      vec res = X2b - 2*X*(b%X_avg) + as_scalar(sum(X_avg % X_avg % b));
      return(res);
    }
    
    vec compute_X2ty(const vec& y, const vec& X_avg) {
      vec X2ty(X.n_cols); X2ty = trans((y.t() * X2));
      vec res = X2ty - 2*trans((y.t() * X))%X_avg + ((X_avg % X_avg) * sum(y));
      return(res);
    }
  
  protected:
    // I think SparseMatrixBlock can overload these with sp_mat types?
    mat X;
    mat X2;
};

// Dense Numeric matrix block
class DenseNumericMatrixBlock: public NumericMatrixBlock {
  public:
    DenseNumericMatrixBlock(mat X_in) {
      X = X_in;
      X2 = X%X;
      nrow = X.n_rows;
      ncol = X.n_cols;
    }
    
    ~DenseNumericMatrixBlock() {}
    
  // protected:
  //   mat X; // data
  //   mat X2; // X%X
};

// Sparse Numeric matrix block
class SparseNumericMatrixBlock: public NumericMatrixBlock {
  public:
    SparseNumericMatrixBlock(sp_mat X_in) {
      X = X_in;
      X2 = X%X;
      nrow = X.n_rows;
      ncol = X.n_cols;
    }
    
    ~SparseNumericMatrixBlock() {}
    
  protected:
    sp_mat X; // data
    sp_mat X2; // X%X
};

// function to take in a sorted vec cuts and tells you which bucket you'd fall it (right-closed)
// faster, if we already have the sort order of t, which we need anyway
array<uvec, 2> get_buckets_and_counts(const vec t, const uvec order_t, const vec cuts) {
  uvec counts(cuts.size()+1, fill::zeros);
  uvec t_to_bin(t.size());
  unsigned int j = 0;
  for (unsigned int i = 0; i < cuts.size(); i++) {
    while (t[order_t[j]] <= cuts[i]) {
      counts[i]++;
      t_to_bin[order_t[j]] = i;
      j++;
    }
  }
  unsigned int i = cuts.size();
  while (j < t_to_bin.size()) {
    counts[i]++;
    t_to_bin[order_t[j]] = i;
    j++;
  }
  
  array<uvec, 2> res;
  res[0] = t_to_bin; res[1] = counts;
  return(res);
}

// faster version of reverse(cumsum(reverse(x)))
vec rev_cumsum_rev(vec x) {
  vec res(x.size());
  res[res.size()-1] = x[x.size()-1];
  for (size_t i=res.size()-1; i>0; i--) {
    res[i-1] = res[i] + x[i-1];
  }
  return(res);
}


// Trend filtering matrix block, i.e. 'stump' terms
class TFGMatrixBlock: public MatrixBlock {
  public:
    TFGMatrixBlock(const vec t, const vec br_in) {
      br = unique(br_in); // gets unique, and sorts for us
      ncol = br.n_rows + 1;
      nrow = t.n_rows;
      
      uvec order_t = sort_index(t);
      
      array<uvec, 2> buckets_and_counts = get_buckets_and_counts(t, order_t, br);
      
      null_bin = index_max(buckets_and_counts[1]);
      bin_to_t = cumsum(buckets_and_counts[1]);
      uvec nz = find(buckets_and_counts[0] != null_bin);
      if (2*nz.size() < nrow) { // if more efficient to store sparse, do so
        t_to_bin = join_cols(buckets_and_counts[0].elem(nz), nz);
      } else { // otherwise, store normally
        t_to_bin = umat(buckets_and_counts[0]);
      }
      
      // now deal with ordering below and above null bin
      // if ((null_bin == 0) || ((null_bin == 1) && (bin_to_t[0] == bin_to_t[1]))) { // corner case
      if ((null_bin == 0) || (bin_to_t[null_bin-1] == 0)) { // corner case
        order_t_low.reset();
      } else {
        order_t_low = order_t.rows(0, bin_to_t[null_bin-1]-1);
      }
      // if ((null_bin == ncol-1) || ((null_bin == ncol-2) && (bin_to_t[bin_to_t.size()-2] == bin_to_t[bin_to_t.size()-1]))) { // corner case
      if (bin_to_t[null_bin] == as_scalar(bin_to_t.tail(1))) { // corner case
        order_t_high.reset();
      } else {
        order_t_high = order_t.rows(bin_to_t[null_bin], order_t.size()-1);
      }
    }
    
    // overload if t is a column from a sparse matrix
    TFGMatrixBlock(const sp_mat t, const vec br_in) {
      br = unique(br_in); // gets unique, and sorts for us
      ncol = br.n_rows + 1;
      nrow = t.n_rows;
      
      // have to make it dense to get sort index.... (MAYBE FIX LATER)
      vec t_dense(t.col(0));
      uvec order_t = sort_index(t_dense); 
      
      array<uvec, 2> buckets_and_counts = get_buckets_and_counts(t_dense, order_t, br);
      
      null_bin = index_max(buckets_and_counts[1]);
      bin_to_t = cumsum(buckets_and_counts[1]);
      uvec nz = find(buckets_and_counts[0] != null_bin);
      if (2*nz.size() < nrow) { // if more efficient to store sparse, do so
        t_to_bin = join_cols(buckets_and_counts[0].elem(nz), nz);
      } else { // otherwise, store normally
        t_to_bin = umat(buckets_and_counts[0]);
      }
      
      // now deal with ordering below and above null bin
      // if ((null_bin == 0) || ((null_bin == 1) && (bin_to_t[0] == bin_to_t[1]))) { // corner case
      if ((null_bin == 0) || (bin_to_t[null_bin-1] == 0)) { // corner case
        order_t_low.reset();
      } else {
        order_t_low = order_t.rows(0, bin_to_t[null_bin-1]-1);
      }
      // if ((null_bin == ncol-1) || ((null_bin == ncol-2) && (bin_to_t[bin_to_t.size()-2] == bin_to_t[bin_to_t.size()-1]))) { // corner case
      if (bin_to_t[null_bin] == as_scalar(bin_to_t.tail(1))) { // corner case
        order_t_high.reset();
      } else {
        order_t_high = order_t.rows(bin_to_t[null_bin], order_t.size()-1);
      }
    }
    
    ~TFGMatrixBlock() {}
    
    vec compute_Xb(const vec& b, const vec& X_avg) {
      vec rev_csb = rev_cumsum_rev(b);
      vec res(nrow);
      if (t_to_bin.n_cols == 1) {
        res = rev_csb.elem(t_to_bin.col(0));
      } else {
        res.fill(rev_csb[null_bin]);
        res.elem(t_to_bin.col(1)) = rev_csb.elem(t_to_bin.col(0));
      }
      res = res - as_scalar(dot(b, X_avg));
      return(res);
    }
    
    vec compute_Xty(const vec& y, const vec& X_avg) {
      vec csy(nrow);
      if (t_to_bin.n_cols == 1) {
        csy = cumsum(y.elem(join_cols(order_t_low, find(t_to_bin.col(0) == null_bin), order_t_high)));
      } else {
        uvec one_to_n = regspace<uvec>(0, nrow-1);
        vector<int> diff_i;
        set_difference(one_to_n.begin(), one_to_n.end(), t_to_bin.col(1).begin(), t_to_bin.col(1).end(),
                       inserter(diff_i, diff_i.begin()));
        uvec which_null_bin = conv_to<uvec>::from(diff_i);
        csy = cumsum(y.elem(join_cols(order_t_low, which_null_bin, order_t_high)));
      }
      vec res = csy.elem(bin_to_t - 1);
      if (bin_to_t[0] == 0) { // weird corner-case in R, CHECK IF IT'S NEEDED IN C++
        vec z = {0.0};
        res = join_cols(z, res);
      }
      res = res - (X_avg  * sum(y));
      return(res);
    }
    
    vec compute_X2b(const vec& b, const vec& X_avg) {
      vec b2 = b%(1 - 2*X_avg);
      vec z(X_avg.n_rows, fill::zeros);
      vec res = compute_Xb(b2, z) + as_scalar(sum(X_avg % X_avg % b));
      return(res);
    }
    
    vec compute_X2ty(const vec& y, const vec& X_avg) {
      vec z(X_avg.n_rows, fill::zeros);
      vec res = compute_Xty(y, z)%(1 - 2*X_avg) + ((X_avg % X_avg) * sum(y));
      return(res);
    }
  
  protected:
    vec br;
    unsigned int null_bin;
    uvec bin_to_t;
    umat t_to_bin; // either just a vector, or a matrix with 2 columns (first is x-values of non-null bin, second is index of non-null bin values)
    uvec order_t_low = {0};
    uvec order_t_high = {0};
};

class StumpsMatrix {
  public:
    // constructor when using dense X
    StumpsMatrix(mat X, uvec _include_linear, uvec _include_stumps, vector<vec> cuts , unsigned int _ncores) {
      ncores = _ncores;
      include_linear = _include_linear;
      include_stumps = _include_stumps;
      uvec which_incl_stumps = find(_include_stumps);
      ncol_lin = sum(include_linear);
      //blocks = vector<MatrixBlock> (which_incl_stumps.size(), MatrixBlock());
      bool any_incl_lin = (ncol_lin > 0);
      if (any_incl_lin) {
        blocks = vector< unique_ptr< MatrixBlock > > (which_incl_stumps.size() + 1);
      } else {
        blocks = vector< unique_ptr< MatrixBlock > > (which_incl_stumps.size());
      }
      
      // first, make stumps blocks (in parallel)
      #if defined(_OPENMP)
        #pragma omp parallel for schedule(dynamic) num_threads(ncores)
      #endif
      for (size_t i = 0; i < which_incl_stumps.size(); i++) {
        // blocks[i] = TFGMatrixBlock(X.col(which_incl_stumps[i]), cuts[which_incl_stumps[i]]);
        blocks[i + any_incl_lin].reset( new TFGMatrixBlock(X.col(which_incl_stumps[i]), cuts[which_incl_stumps[i]]) );
      }
      
      // now, add linear part
      // ADD CODE HERE TO GET COLUMN SCALE FACTORS AND SCALE X (need extra input to constructor for, e.g. 'sd', or 'max' or 'none' scaling)
      if (any_incl_lin) {
        blocks.front().reset( new DenseNumericMatrixBlock(X.cols(find(include_linear))) );
      }
      
      nrow = X.n_rows;
      ncol = 0;
      // get ncol, and col_offsets for each block
      for (size_t i = 0; i < blocks.size(); i++) {
        blocks[i].get()->col_offset = ncol;
        ncol += blocks[i].get()->ncol;
      }
    }
    
    // overload constructor when using sparse X
    StumpsMatrix(sp_mat X, uvec _include_linear, uvec _include_stumps, vector<vec> cuts , unsigned int _ncores) {
      ncores = _ncores;
      include_linear = _include_linear;
      include_stumps = _include_stumps;
      ncol_lin = sum(include_linear);
      uvec which_incl_stumps = find(_include_stumps);
      //blocks = vector<MatrixBlock> (which_incl_stumps.size(), MatrixBlock());
      bool any_incl_lin = (ncol_lin > 0);
      if (any_incl_lin) {
        blocks = vector< unique_ptr< MatrixBlock > > (which_incl_stumps.size() + 1);
      } else {
        blocks = vector< unique_ptr< MatrixBlock > > (which_incl_stumps.size());
      }
      
      // first, make stumps blocks (in parallel)
      #if defined(_OPENMP)
        #pragma omp parallel for schedule(dynamic) num_threads(ncores)
      #endif
      for (size_t i = 0; i < which_incl_stumps.size(); i++) {
        // blocks[i] = TFGMatrixBlock(X.col(which_incl_stumps[i]), cuts[which_incl_stumps[i]]);
        blocks[i + any_incl_lin].reset( new TFGMatrixBlock(X.col(which_incl_stumps[i]), cuts[which_incl_stumps[i]]) );
      }
      
      // now, add linear part
      if (any_incl_lin) {
        blocks.front().reset( new SparseNumericMatrixBlock(X.cols(find(include_linear))) );
      }
      
      nrow = X.n_rows;
      ncol = 0;
      // get ncol, and col_offsets for each block
      for (size_t i = 0; i < blocks.size(); i++) {
        blocks[i].get()->col_offset = ncol;
        ncol += blocks[i].get()->ncol;
      }
    }
    
    ~StumpsMatrix() {
      for (size_t i = 0; i < blocks.size(); i++) {
        blocks[i].reset();
      }
    }

    #if defined(_OPENMP)
      #pragma omp declare reduction( + : vec : omp_out += omp_in ) \
        initializer( omp_priv = zeros<vec>(omp_orig.n_rows))
    #endif
    
    vec compute_Xb(const vec& b, const vec& X_avg) {
      vec Xb(nrow, fill::zeros);
      #if defined(_OPENMP)
        #pragma omp parallel for reduction(+:Xb) schedule(dynamic) num_threads(ncores)
      #endif
      for(size_t i = 0; i < blocks.size(); i++) {
        vec my_Xb = blocks[i].get()->compute_Xb(b.rows(blocks[i].get()->col_offset, blocks[i].get()->col_offset + blocks[i].get()->ncol - 1), X_avg.rows(blocks[i].get()->col_offset, blocks[i].get()->col_offset + blocks[i].get()->ncol - 1));
        Xb += my_Xb;
      }
      return(Xb);
    }
    
    vec compute_Xty(const vec& y, const vec& X_avg) {
      vec Xty(ncol);
      #if defined(_OPENMP)
        #pragma omp parallel for schedule(dynamic) num_threads(ncores)
      #endif
      for(size_t i = 0; i < blocks.size(); i++) {
        Xty.rows(blocks[i].get()->col_offset, blocks[i].get()->col_offset + blocks[i].get()->ncol - 1) = blocks[i].get()->compute_Xty(y, X_avg.rows(blocks[i].get()->col_offset, blocks[i].get()->col_offset + blocks[i].get()->ncol - 1));
      }
      return(Xty);
    }
    
    vec compute_X2b(const vec& b, const vec& X_avg) {
      vec X2b(nrow, fill::zeros);
      #if defined(_OPENMP)
        #pragma omp parallel for reduction(+:X2b) schedule(dynamic) num_threads(ncores)
      #endif
      for(size_t i = 0; i < blocks.size(); i++) {
        vec my_X2b = blocks[i].get()->compute_X2b(b.rows(blocks[i].get()->col_offset, blocks[i].get()->col_offset + blocks[i].get()->ncol - 1), X_avg.rows(blocks[i].get()->col_offset, blocks[i].get()->col_offset + blocks[i].get()->ncol - 1));
        X2b += my_X2b;
      }
      return(X2b);
    }
    
    vec compute_X2ty(const vec& y, const vec& X_avg) {
      vec X2ty(ncol);
      #if defined(_OPENMP)
        #pragma omp parallel for schedule(dynamic) num_threads(ncores)
      #endif
      for(size_t i = 0; i < blocks.size(); i++) {
        X2ty.rows(blocks[i].get()->col_offset, blocks[i].get()->col_offset + blocks[i].get()->ncol - 1) = blocks[i].get()->compute_X2ty(y, X_avg.rows(blocks[i].get()->col_offset, blocks[i].get()->col_offset + blocks[i].get()->ncol - 1));
      }
      return(X2ty);
    }
    
    unsigned int nrow = 0;
    unsigned int ncol = 0;
    uvec include_linear;
    uvec include_stumps;
    vec col_scale_factors;
    unsigned int ncol_lin;
    
  protected:
    vector< unique_ptr< MatrixBlock > > blocks;
    unsigned int ncores = 1;
};

RCPP_MODULE(stumpsMatrix_module) {
  
  class_<StumpsMatrix>( "StumpsMatrix" )
  
  .constructor<mat, uvec, uvec, vector<vec> , unsigned int>()
  .constructor<sp_mat, uvec, uvec, vector<vec> , unsigned int>()
  
  .field( "nrow", &StumpsMatrix::nrow )
  .field( "ncol", &StumpsMatrix::ncol )
  .field( "include_linear", &StumpsMatrix::include_linear )
  .field( "include_stumps", &StumpsMatrix::include_stumps )
  .field( "col_scale_factors", &StumpsMatrix::col_scale_factors )
}


#endif
