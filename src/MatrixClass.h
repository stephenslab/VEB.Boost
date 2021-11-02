#ifndef INCLUDE_MatrixClass_h
#define INCLUDE_MatrixClass_h

#include <RcppArmadillo.h>
#include <memory>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

namespace stumpsmatrix {
  // Matrix block abstract class, i.e. an element of X_stumps (either numeric matrix, or tfg matrix)
  class MatrixBlock {
    public: 
      MatrixBlock() {}
      
      ~MatrixBlock() {}
      
      // method to compute (X - X_avg) * b
      virtual arma::vec compute_Xb(const arma::vec& b, const arma::vec& X_avg) = 0;
      
      // method to compute (X - X_avg)' * y
      virtual arma::vec compute_Xty(const arma::vec& y, const arma::vec& X_avg) = 0;
      
      // method to compute (X - X_avg)^2 * b
      virtual arma::vec compute_X2b(const arma::vec& b, const arma::vec& X_avg) = 0;
      
      // method to compute (X' - X_avg')^2 * y
      virtual arma::vec compute_X2ty(const arma::vec& y, const arma::vec& X_avg) = 0;
      
      unsigned int nrow, ncol;
      unsigned int col_offset; // column offset in entire matrix, i.e. sum of # columns before this block
  };
  
  // Numeric matrix block, i.e. first element of X_stumps when including linear
  class NumericMatrixBlock: public MatrixBlock {
    public:
      ~NumericMatrixBlock() {}
      
      arma::vec compute_Xb(const arma::vec& b, const arma::vec& X_avg) {
        arma::vec res = X*b - arma::as_scalar(arma::dot(b, X_avg));
        return(res);
      }
      
      arma::vec compute_Xty(const arma::vec& y, const arma::vec& X_avg) {
        arma::vec res = arma::trans((y.t() * X)) - (X_avg  * arma::sum(y));
        return(res);
      }
      
      arma::vec compute_X2b(const arma::vec& b, const arma::vec& X_avg) {
        arma::vec X2b(X.n_rows); X2b = X2*b;
        arma::vec res = X2b - 2*X*(b%X_avg) + arma::as_scalar(arma::sum(X_avg % X_avg % b));
        return(res);
      }
      
      arma::vec compute_X2ty(const arma::vec& y, const arma::vec& X_avg) {
        arma::vec X2ty(X.n_cols); X2ty = arma::trans((y.t() * X2));
        arma::vec res = X2ty - 2*arma::trans((y.t() * X))%X_avg + ((X_avg % X_avg) * arma::sum(y));
        return(res);
      }
    
    protected:
      // I think SparseMatrixBlock can overload these with arma::sp_mat types?
      arma::mat X;
      arma::mat X2;
  };
  
  // Dense Numeric matrix block
  class DenseNumericMatrixBlock: public NumericMatrixBlock {
    public:
      DenseNumericMatrixBlock(arma::mat X_in) {
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
      SparseNumericMatrixBlock(arma::sp_mat X_in) {
        X = X_in;
        X2 = X%X;
        nrow = X.n_rows;
        ncol = X.n_cols;
      }
      
      ~SparseNumericMatrixBlock() {}
      
    protected:
      arma::sp_mat X; // data
      arma::sp_mat X2; // X%X
  };
  
  // function to take in a sorted vec cuts and tells you which bucket you'd fall it (right-closed)
  // faster, if we already have the sort order of t, which we need anyway
  std::array<arma::uvec, 2> get_buckets_and_counts(const arma::vec& t, const arma::uvec& order_t, const arma::vec& cuts) {
    arma::uvec counts(cuts.size()+1, arma::fill::zeros);
    arma::uvec t_to_bin(t.size());
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
    
    std::array<arma::uvec, 2> res;
    res[0] = t_to_bin; res[1] = counts;
    return(res);
  }
  
  // faster version of reverse(cumsum(reverse(x)))
  arma::vec rev_cumsum_rev(const arma::vec& x) {
    arma::vec res(x.size());
    res[res.size()-1] = x[x.size()-1];
    for (size_t i=res.size()-1; i>0; i--) {
      res[i-1] = res[i] + x[i-1];
    }
    return(res);
  }
  
  
  // Trend filtering matrix block, i.e. 'stump' terms
  class TFGMatrixBlock: public MatrixBlock {
    public:
      TFGMatrixBlock(const arma::vec& t, const arma::vec& br_in) {
        br = arma::unique(br_in); // gets unique, and sorts for us
        ncol = br.n_rows + 1;
        nrow = t.n_rows;
        
        arma::uvec order_t = arma::sort_index(t);
        
        std::array<arma::uvec, 2> buckets_and_counts = get_buckets_and_counts(t, order_t, br);
        
        null_bin = arma::index_max(buckets_and_counts[1]);
        bin_to_t = arma::cumsum(buckets_and_counts[1]);
        arma::uvec nz = arma::find(buckets_and_counts[0] != null_bin);
        if (2*nz.size() < nrow) { // if more efficient to store sparse, do so
          t_to_bin = arma::join_cols(buckets_and_counts[0].elem(nz), nz);
        } else { // otherwise, store normally
          t_to_bin = arma::umat(buckets_and_counts[0]);
        }
        
        // now deal with ordering below and above null bin
        // if ((null_bin == 0) || ((null_bin == 1) && (bin_to_t[0] == bin_to_t[1]))) { // corner case
        if ((null_bin == 0) || (bin_to_t[null_bin-1] == 0)) { // corner case
          order_t_low.reset();
        } else {
          order_t_low = order_t.rows(0, bin_to_t[null_bin-1]-1);
        }
        // if ((null_bin == ncol-1) || ((null_bin == ncol-2) && (bin_to_t[bin_to_t.size()-2] == bin_to_t[bin_to_t.size()-1]))) { // corner case
        if (bin_to_t[null_bin] == arma::as_scalar(bin_to_t.tail(1))) { // corner case
          order_t_high.reset();
        } else {
          order_t_high = order_t.rows(bin_to_t[null_bin], order_t.size()-1);
        }
      }
      
      // overload if t is a column from a sparse matrix
      TFGMatrixBlock(const arma::sp_mat& t, const arma::vec& br_in) {
        br = arma::unique(br_in); // gets unique, and sorts for us
        ncol = br.n_rows + 1;
        nrow = t.n_rows;
        
        // have to make it dense to get sort index.... (MAYBE FIX LATER)
        arma::vec t_dense(t.col(0));
        arma::uvec order_t = arma::sort_index(t_dense); 
        
        std::array<arma::uvec, 2> buckets_and_counts = get_buckets_and_counts(t_dense, order_t, br);
        
        null_bin = arma::index_max(buckets_and_counts[1]);
        bin_to_t = arma::cumsum(buckets_and_counts[1]);
        arma::uvec nz = arma::find(buckets_and_counts[0] != null_bin);
        if (2*nz.size() < nrow) { // if more efficient to store sparse, do so
          t_to_bin = arma::join_cols(buckets_and_counts[0].elem(nz), nz);
        } else { // otherwise, store normally
          t_to_bin = arma::umat(buckets_and_counts[0]);
        }
        
        // now deal with ordering below and above null bin
        // if ((null_bin == 0) || ((null_bin == 1) && (bin_to_t[0] == bin_to_t[1]))) { // corner case
        if ((null_bin == 0) || (bin_to_t[null_bin-1] == 0)) { // corner case
          order_t_low.reset();
        } else {
          order_t_low = order_t.rows(0, bin_to_t[null_bin-1]-1);
        }
        // if ((null_bin == ncol-1) || ((null_bin == ncol-2) && (bin_to_t[bin_to_t.size()-2] == bin_to_t[bin_to_t.size()-1]))) { // corner case
        if (bin_to_t[null_bin] == arma::as_scalar(bin_to_t.tail(1))) { // corner case
          order_t_high.reset();
        } else {
          order_t_high = order_t.rows(bin_to_t[null_bin], order_t.size()-1);
        }
      }
      
      ~TFGMatrixBlock() {}
      
      arma::vec compute_Xb(const arma::vec& b, const arma::vec& X_avg) {
        arma::vec rev_csb = rev_cumsum_rev(b);
        arma::vec res(nrow);
        if (t_to_bin.n_cols == 1) {
          res = rev_csb.elem(t_to_bin.col(0));
        } else {
          res.fill(rev_csb[null_bin]);
          res.elem(t_to_bin.col(1)) = rev_csb.elem(t_to_bin.col(0));
        }
        res = res - arma::as_scalar(dot(b, X_avg));
        return(res);
      }
      
      arma::vec compute_Xty(const arma::vec& y, const arma::vec& X_avg) {
        arma::vec csy(nrow);
        if (t_to_bin.n_cols == 1) {
          csy = arma::cumsum(y.elem(arma::join_cols(order_t_low, arma::find(t_to_bin.col(0) == null_bin), order_t_high)));
        } else {
          arma::uvec one_to_n = arma::regspace<arma::uvec>(0, nrow-1);
          std::vector<int> diff_i;
          std::set_difference(one_to_n.begin(), one_to_n.end(), t_to_bin.col(1).begin(), t_to_bin.col(1).end(),
                         std::inserter(diff_i, diff_i.begin()));
          arma::uvec which_null_bin = arma::conv_to<arma::uvec>::from(diff_i);
          csy = arma::cumsum(y.elem(arma::join_cols(order_t_low, which_null_bin, order_t_high)));
        }
        arma::vec res = csy.elem(bin_to_t - 1);
        if (bin_to_t[0] == 0) { // weird corner-case in R, CHECK IF IT'S NEEDED IN C++
          arma::vec z = {0.0};
          res = arma::join_cols(z, res);
        }
        res = res - (X_avg  * arma::sum(y));
        return(res);
      }
      
      arma::vec compute_X2b(const arma::vec& b, const arma::vec& X_avg) {
        arma::vec b2 = b%(1 - 2*X_avg);
        arma::vec z(X_avg.n_rows, arma::fill::zeros);
        arma::vec res = compute_Xb(b2, z) + arma::as_scalar(arma::sum(X_avg % X_avg % b));
        return(res);
      }
      
      arma::vec compute_X2ty(const arma::vec& y, const arma::vec& X_avg) {
        arma::vec z(X_avg.n_rows, arma::fill::zeros);
        arma::vec res = compute_Xty(y, z)%(1 - 2*X_avg) + ((X_avg % X_avg) * arma::sum(y));
        return(res);
      }
    
    protected:
      arma::vec br;
      unsigned int null_bin;
      arma::uvec bin_to_t;
      arma::umat t_to_bin; // either just a std::vector, or a matrix with 2 columns (first is x-values of non-null bin, second is index of non-null bin values)
      arma::uvec order_t_low = {0};
      arma::uvec order_t_high = {0};
  };
  
  class StumpsMatrix {
    public:
      // constructor when using dense X
      StumpsMatrix(arma::mat X, arma::uvec _include_linear, arma::uvec _include_stumps, std::vector<arma::vec> cuts , unsigned int _ncores) {
        ncores = _ncores;
        include_linear = _include_linear;
        include_stumps = _include_stumps;
        arma::uvec which_incl_stumps = arma::find(_include_stumps);
        ncol_lin = arma::sum(include_linear);
        //blocks = std::vector<MatrixBlock> (which_incl_stumps.size(), MatrixBlock());
        bool any_incl_lin = (ncol_lin > 0);
        if (any_incl_lin) {
          blocks = std::vector< std::unique_ptr< MatrixBlock > > (which_incl_stumps.size() + 1);
        } else {
          blocks = std::vector< std::unique_ptr< MatrixBlock > > (which_incl_stumps.size());
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
          blocks.front().reset( new DenseNumericMatrixBlock(X.cols(arma::find(include_linear))) );
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
      StumpsMatrix(arma::sp_mat X, arma::uvec _include_linear, arma::uvec _include_stumps, std::vector<arma::vec> cuts , unsigned int _ncores) {
        ncores = _ncores;
        include_linear = _include_linear;
        include_stumps = _include_stumps;
        ncol_lin = sum(include_linear);
        arma::uvec which_incl_stumps = arma::find(_include_stumps);
        //blocks = std::vector<MatrixBlock> (which_incl_stumps.size(), MatrixBlock());
        bool any_incl_lin = (ncol_lin > 0);
        if (any_incl_lin) {
          blocks = std::vector< std::unique_ptr< MatrixBlock > > (which_incl_stumps.size() + 1);
        } else {
          blocks = std::vector< std::unique_ptr< MatrixBlock > > (which_incl_stumps.size());
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
          blocks.front().reset( new SparseNumericMatrixBlock(X.cols(arma::find(include_linear))) );
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
        #pragma omp declare reduction( + : arma::vec : omp_out += omp_in ) \
          initializer( omp_priv = arma::zeros<arma::vec>(omp_orig.n_rows))
      #endif
      
      arma::vec compute_Xb(const arma::vec& b, const arma::vec& X_avg) {
        arma::vec Xb(nrow, arma::fill::zeros);
        #if defined(_OPENMP)
          #pragma omp parallel for reduction(+:Xb) schedule(dynamic) num_threads(ncores)
        #endif
        for(size_t i = 0; i < blocks.size(); i++) {
          arma::vec my_Xb = blocks[i].get()->compute_Xb(b.rows(blocks[i].get()->col_offset, blocks[i].get()->col_offset + blocks[i].get()->ncol - 1), X_avg.rows(blocks[i].get()->col_offset, blocks[i].get()->col_offset + blocks[i].get()->ncol - 1));
          Xb += my_Xb;
        }
        return(Xb);
      }
      
      arma::vec compute_Xty(const arma::vec& y, const arma::vec& X_avg) {
        arma::vec Xty(ncol);
        #if defined(_OPENMP)
          #pragma omp parallel for schedule(dynamic) num_threads(ncores)
        #endif
        for(size_t i = 0; i < blocks.size(); i++) {
          Xty.rows(blocks[i].get()->col_offset, blocks[i].get()->col_offset + blocks[i].get()->ncol - 1) = blocks[i].get()->compute_Xty(y, X_avg.rows(blocks[i].get()->col_offset, blocks[i].get()->col_offset + blocks[i].get()->ncol - 1));
        }
        return(Xty);
      }
      
      arma::vec compute_X2b(const arma::vec& b, const arma::vec& X_avg) {
        arma::vec X2b(nrow, arma::fill::zeros);
        #if defined(_OPENMP)
          #pragma omp parallel for reduction(+:X2b) schedule(dynamic) num_threads(ncores)
        #endif
        for(size_t i = 0; i < blocks.size(); i++) {
          arma::vec my_X2b = blocks[i].get()->compute_X2b(b.rows(blocks[i].get()->col_offset, blocks[i].get()->col_offset + blocks[i].get()->ncol - 1), X_avg.rows(blocks[i].get()->col_offset, blocks[i].get()->col_offset + blocks[i].get()->ncol - 1));
          X2b += my_X2b;
        }
        return(X2b);
      }
      
      arma::vec compute_X2ty(const arma::vec& y, const arma::vec& X_avg) {
        arma::vec X2ty(ncol);
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
      arma::uvec include_linear;
      arma::uvec include_stumps;
      arma::vec col_scale_factors;
      unsigned int ncol_lin;
      
    protected:
      std::vector< std::unique_ptr< MatrixBlock > > blocks;
      unsigned int ncores = 1;
  };
  
}

// RCPP_MODULE(stumpsMatrix_module) {
//   
//   class_<StumpsMatrix>( "StumpsMatrix" )
//   
//   .constructor<arma::mat, arma::uvec, arma::uvec, std::vector<arma::vec> , unsigned int>()
//   .constructor<arma::sp_mat, arma::uvec, arma::uvec, std::vector<arma::vec> , unsigned int>()
//   
//   .field( "nrow", &StumpsMatrix::nrow )
//   .field( "ncol", &StumpsMatrix::ncol )
//   .field( "include_linear", &StumpsMatrix::include_linear )
//   .field( "include_stumps", &StumpsMatrix::include_stumps )
//   .field( "col_scale_factors", &StumpsMatrix::col_scale_factors )
// }


#endif
