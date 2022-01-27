#pragma once

#include <RcppArmadillo.h>
#include <memory>
#include <array>
#if defined(_OPENMP)
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

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

    //protected:
      arma::vec col_scale_factors; // 1 / column scaling factors (so multiply columns by these factors)
      arma::uvec include_linear;
      arma::vec br;
      unsigned int null_bin;
      arma::uvec bin_to_t;
      arma::umat t_to_bin; // either just a std::vector, or a matrix with 2 columns (first is x-values of non-null bin, second is index of non-null bin values)
      arma::uvec order_t_low = { 0 };
      arma::uvec order_t_high = { 0 };
  };
  
  //// Numeric matrix block, i.e. first element of X_stumps when including linear
  //class NumericMatrixBlock: public MatrixBlock {
  //  public:
  //    ~NumericMatrixBlock() {}
  //    
  //    arma::vec compute_Xb(const arma::vec& b, const arma::vec& X_avg) {
  //      arma::vec res = X*b - arma::as_scalar(arma::dot(b, X_avg));
  //      return(res);
  //    }
  //    
  //    arma::vec compute_Xty(const arma::vec& y, const arma::vec& X_avg) {
  //      arma::vec res = arma::trans((y.t() * X)) - (X_avg  * arma::sum(y));
  //      return(res);
  //    }
  //    
  //    arma::vec compute_X2b(const arma::vec& b, const arma::vec& X_avg) {
  //      arma::vec X2b(X.n_rows); X2b = X2*b;
  //      arma::vec res = X2b - 2*X*(b%X_avg) + arma::as_scalar(arma::sum(X_avg % X_avg % b));
  //      return(res);
  //    }
  //    
  //    arma::vec compute_X2ty(const arma::vec& y, const arma::vec& X_avg) {
  //      arma::vec X2ty(X.n_cols); X2ty = arma::trans((y.t() * X2));
  //      arma::vec res = X2ty - 2*arma::trans((y.t() * X))%X_avg + ((X_avg % X_avg) * arma::sum(y));
  //      return(res);
  //    }
  //  
  //  protected:
  //    // I think SparseMatrixBlock can overload these with arma::sp_mat types?
  //    arma::mat X;
  //    arma::mat X2;
  //};
  
  // Dense Numeric matrix block
  class DenseNumericMatrixBlock: public MatrixBlock {
    public:
      DenseNumericMatrixBlock(arma::mat& X_in, arma::uvec include_linear_in, arma::vec col_scale_factors_in) {
        X = arma::mat(X_in.memptr(), X_in.n_rows, X_in.n_cols, false, false);
        //X = X_in;
        include_linear = arma::find(include_linear_in);
        col_scale_factors = 1 / col_scale_factors_in;
        X2 = arma::square(X.cols(include_linear) * arma::diagmat(col_scale_factors));
        nrow = X.n_rows;
        ncol = X.n_cols;
      }

      // input is existing DenseNumericMatrixBlock
      DenseNumericMatrixBlock(arma::mat& X_in, std::shared_ptr<MatrixBlock> Mat_in) {
        X = arma::mat(X_in.memptr(), X_in.n_rows, X_in.n_cols, false, false);
        include_linear = Mat_in.get()->include_linear;
        col_scale_factors = Mat_in.get()->col_scale_factors;
        X2 = arma::square(X.cols(include_linear) * arma::diagmat(col_scale_factors));
        nrow = X.n_rows;
        ncol = X.n_cols;
      }
      
      ~DenseNumericMatrixBlock() {}

      arma::vec compute_Xb(const arma::vec& b, const arma::vec& X_avg) {
          arma::vec res = X.cols(include_linear) * (b % col_scale_factors) - arma::as_scalar(arma::dot(b, X_avg));
          return(res);
      }

      arma::vec compute_Xty(const arma::vec& y, const arma::vec& X_avg) {
          arma::vec res = (arma::trans(y.t() * X.cols(include_linear)) % col_scale_factors) - (X_avg * arma::sum(y));
          return(res);
      }

      arma::vec compute_X2b(const arma::vec& b, const arma::vec& X_avg) {
          arma::vec res = X2 * b - 2 * X.cols(include_linear) * (b % X_avg % col_scale_factors) + arma::as_scalar(arma::sum(X_avg % X_avg % b));
          return(res);
      }

      arma::vec compute_X2ty(const arma::vec& y, const arma::vec& X_avg) {
          arma::vec res = arma::trans(y.t() * X2) - 2 * arma::trans(y.t() * X.cols(include_linear)) % (X_avg % col_scale_factors) + ((X_avg % X_avg) * arma::sum(y));
          return(res);
      }
      
    protected:
        arma::mat X; // data
        arma::mat X2; // scaled_data % scaled_data (have to compute when multiplying anyway, so might as well store it)
        //arma::vec col_scale_factors; // 1 / column scaling factors (so multiply columns by these factors)
        //arma::uvec include_linear;
  };
  
  // Sparse Numeric matrix block
  class SparseNumericMatrixBlock: public MatrixBlock {
    public:
      SparseNumericMatrixBlock(arma::sp_mat& X_in, arma::uvec include_linear_in, arma::vec col_scale_factors_in) {
        X = X_in;
        include_linear = arma::find(include_linear_in);
        col_scale_factors = 1 / col_scale_factors_in;
        X2 = arma::square(X.cols(include_linear) * arma::diagmat(col_scale_factors));
        nrow = X.n_rows;
        ncol = X.n_cols;
      }

      // input is existing SparseNumericMatrixBlock
      SparseNumericMatrixBlock(arma::sp_mat& X_in, std::shared_ptr<MatrixBlock> Mat_in) {
          X = X_in;
          include_linear = Mat_in.get()->include_linear;
          col_scale_factors = Mat_in.get()->col_scale_factors;
          X2 = arma::square(X.cols(include_linear) * arma::diagmat(col_scale_factors));
          nrow = X.n_rows;
          ncol = X.n_cols;
      }
      
      ~SparseNumericMatrixBlock() {}

      arma::vec compute_Xb(const arma::vec& b, const arma::vec& X_avg) {
          arma::vec res = X.cols(include_linear) * (b % col_scale_factors) - arma::as_scalar(arma::dot(b, X_avg));
          return(res);
      }

      arma::vec compute_Xty(const arma::vec& y, const arma::vec& X_avg) {
          arma::vec res = (arma::trans(y.t() * X.cols(include_linear)) % col_scale_factors) - (X_avg * arma::sum(y));
          return(res);
      }

      arma::vec compute_X2b(const arma::vec& b, const arma::vec& X_avg) {
          arma::vec res = X2 * b - 2 * X.cols(include_linear) * (b % X_avg % col_scale_factors) + arma::as_scalar(arma::sum(X_avg % X_avg % b));
          return(res);
      }

      arma::vec compute_X2ty(const arma::vec& y, const arma::vec& X_avg) {
          arma::vec res = arma::trans(y.t() * X2) - 2 * arma::trans(y.t() * X.cols(include_linear)) % (X_avg % col_scale_factors) + ((X_avg % X_avg) * arma::sum(y));
          return(res);
      }
      
    protected:
        arma::mat X; // data
        arma::mat X2; // scaled_data % scaled_data (have to compute when multiplying anyway, so might as well store it)
        //arma::vec col_scale_factors; // 1 / column scaling factors (so multiply columns by these factors)
        //arma::uvec include_linear;
  };
  
  // function to take in a sorted vec cuts and tells you which bucket you'd fall it (right-closed)
  // faster, if we already have the sort order of t, which we need anyway
  inline std::array<arma::uvec, 2> get_buckets_and_counts(const arma::vec& t, const arma::uvec& order_t, const arma::vec& cuts) {
    arma::uvec counts(cuts.size()+1, arma::fill::zeros);
    arma::uvec t_to_bin(t.size());
    unsigned int j = 0;
    for (unsigned int i = 0; i < cuts.size(); i++) {
      while ((j < order_t.size()) && (t[order_t[j]] <= cuts[i])) {
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
  inline arma::vec rev_cumsum_rev(const arma::vec& x) {
    arma::vec res(x.size());
    res[res.size()-1] = x[x.size()-1];
    for (size_t i=res.size()-1; i>0; i--) {
      res[i-1] = res[i] + x[i-1];
    }
    return(res);
  }

  // C++ version of ppoints
  inline arma::vec ppoints_cpp(unsigned int num_cuts) {
      double a;
      if (num_cuts <= 10) {
          a = 0.375;
      } else {
          a = 0.5;
      }
      return (arma::linspace(1, num_cuts, num_cuts) - a) / (num_cuts + 1 - 2 * a);
  }

  // C++ version of quantile type 7 (R does type 7, armadillo does type 5)
  
  inline arma::vec quantile_type7(const arma::vec& t, const arma::vec& probs) {
      int n = t.size();
      arma::vec index = std::max(n - 1, 0) * probs;
      arma::uvec lo = arma::conv_to<arma::uvec>::from(arma::floor(index));
      arma::uvec hi = arma::conv_to<arma::uvec>::from(arma::ceil(index));
      arma::vec x = arma::sort(t);
      arma::vec qs = x.elem(lo);
      arma::uvec i = arma::find((index > lo) && (x.elem(hi) != qs));
      arma::vec h = arma::conv_to<arma::vec>::from(index - lo);
      h = h.elem(i);
      qs.elem(i) = (1 - h) % qs.elem(i) + h % x.elem(hi.elem(i));
      return qs;
  }
  

  // functions to get cuts (N.B. for quantiles, armadillo uses R quantile formula 5 instead of 7, also might not work with sp_mat?)
  // get cuts in C++
  inline arma::vec get_cuts_cpp(const arma::vec& t, const unsigned int num_cuts, bool use_quants) {
    arma::vec cuts;
    if (num_cuts == 0) {
        cuts = t;
    } else {
        arma::vec ppts = ppoints_cpp(num_cuts);
        if (use_quants) {
            //cuts = arma::quantile(t, ppts);
            cuts = quantile_type7(t, ppts);
        } else {
            double t_min = arma::min(t);
            double t_max = arma::max(t);
            cuts = t_min * (1 - ppts) + t_max * ppts;
        }
    }
    return cuts;
  }

  // get cuts in C++ (overload for sparse vector)
  inline arma::vec get_cuts_cpp(const arma::sp_mat& t, const unsigned int num_cuts, bool use_quants) {
      arma::vec cuts;
      if (num_cuts == 0) {
          cuts = arma::vec(t.col(0));
      } else {
          arma::vec ppts = ppoints_cpp(num_cuts);
          if (use_quants) {
              // = arma::quantile(arma::vec(t.col(0)), ppts);
              cuts = quantile_type7(arma::vec(t.col(0)), ppts);
          } else {
              double t_min = arma::min(t.col(0));
              double t_max = arma::max(t.col(0));
              cuts = t_min * (1 - ppts) + t_max * ppts;
          }
      }
      return cuts;
  }
  
  // Trend filtering matrix block, i.e. 'stump' terms
  class TFGMatrixBlock: public MatrixBlock {
    public:
      TFGMatrixBlock(const arma::vec& t, unsigned int num_cuts, bool use_quants) {
        br = get_cuts_cpp(t, num_cuts, use_quants);
        br = arma::unique(br); // gets unique, and sorts for us
        ncol = br.n_rows + 1;
        nrow = t.n_rows;
        
        arma::uvec order_t = arma::sort_index(t);
        
        std::array<arma::uvec, 2> buckets_and_counts = get_buckets_and_counts(t, order_t, br);
        
        null_bin = arma::index_max(buckets_and_counts[1]);
        bin_to_t = arma::cumsum(buckets_and_counts[1]);
        arma::uvec nz = arma::find(buckets_and_counts[0] != null_bin);
        if (2*nz.size() < nrow) { // if more efficient to store sparse, do so
          t_to_bin = arma::join_rows(buckets_and_counts[0].elem(nz), nz);
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
      TFGMatrixBlock(const arma::sp_mat& t, unsigned int num_cuts, bool use_quants) {
        br = get_cuts_cpp(t, num_cuts, use_quants);
        br = arma::unique(br); // gets unique, and sorts for us
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
          t_to_bin = arma::join_rows(buckets_and_counts[0].elem(nz), nz);
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

      // input is existing TFGMatrixBlock with dense t
      TFGMatrixBlock(const arma::vec& t, std::shared_ptr<MatrixBlock> TFG_in) {
          br = TFG_in.get()->br;
          ncol = br.n_rows + 1;
          nrow = t.n_rows;

          arma::uvec order_t = arma::sort_index(t);

          std::array<arma::uvec, 2> buckets_and_counts = get_buckets_and_counts(t, order_t, br);

          null_bin = arma::index_max(buckets_and_counts[1]);
          bin_to_t = arma::cumsum(buckets_and_counts[1]);
          arma::uvec nz = arma::find(buckets_and_counts[0] != null_bin);
          if (2 * nz.size() < nrow) { // if more efficient to store sparse, do so
              t_to_bin = arma::join_rows(buckets_and_counts[0].elem(nz), nz);
          }
          else { // otherwise, store normally
              t_to_bin = arma::umat(buckets_and_counts[0]);
          }

          // now deal with ordering below and above null bin
          // if ((null_bin == 0) || ((null_bin == 1) && (bin_to_t[0] == bin_to_t[1]))) { // corner case
          if ((null_bin == 0) || (bin_to_t[null_bin - 1] == 0)) { // corner case
              order_t_low.reset();
          }
          else {
              order_t_low = order_t.rows(0, bin_to_t[null_bin - 1] - 1);
          }
          // if ((null_bin == ncol-1) || ((null_bin == ncol-2) && (bin_to_t[bin_to_t.size()-2] == bin_to_t[bin_to_t.size()-1]))) { // corner case
          if (bin_to_t[null_bin] == arma::as_scalar(bin_to_t.tail(1))) { // corner case
              order_t_high.reset();
          }
          else {
              order_t_high = order_t.rows(bin_to_t[null_bin], order_t.size() - 1);
          }
      }

      // input is existing TFGMatrixBlock with sparse t
      TFGMatrixBlock(const arma::sp_mat& t, std::shared_ptr<MatrixBlock> TFG_in) {
          br = TFG_in.get()->br; // gets unique, and sorts for us
          ncol = br.n_rows + 1;
          nrow = t.n_rows;

          // have to make it dense to get sort index.... (MAYBE FIX LATER)
          arma::vec t_dense(t.col(0));
          arma::uvec order_t = arma::sort_index(t_dense);

          std::array<arma::uvec, 2> buckets_and_counts = get_buckets_and_counts(t_dense, order_t, br);

          null_bin = arma::index_max(buckets_and_counts[1]);
          bin_to_t = arma::cumsum(buckets_and_counts[1]);
          arma::uvec nz = arma::find(buckets_and_counts[0] != null_bin);
          if (2 * nz.size() < nrow) { // if more efficient to store sparse, do so
              t_to_bin = arma::join_rows(buckets_and_counts[0].elem(nz), nz);
          }
          else { // otherwise, store normally
              t_to_bin = arma::umat(buckets_and_counts[0]);
          }

          // now deal with ordering below and above null bin
          // if ((null_bin == 0) || ((null_bin == 1) && (bin_to_t[0] == bin_to_t[1]))) { // corner case
          if ((null_bin == 0) || (bin_to_t[null_bin - 1] == 0)) { // corner case
              order_t_low.reset();
          }
          else {
              order_t_low = order_t.rows(0, bin_to_t[null_bin - 1] - 1);
          }
          // if ((null_bin == ncol-1) || ((null_bin == ncol-2) && (bin_to_t[bin_to_t.size()-2] == bin_to_t[bin_to_t.size()-1]))) { // corner case
          if (bin_to_t[null_bin] == arma::as_scalar(bin_to_t.tail(1))) { // corner case
              order_t_high.reset();
          }
          else {
              order_t_high = order_t.rows(bin_to_t[null_bin], order_t.size() - 1);
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
        //arma::vec res = csy.elem(bin_to_t - 1);
        //if (bin_to_t[0] == 0) { // weird corner-case in R, NEED TO FIX IN C++ (b/c above will index at -1)
        //  arma::vec z = {0.0};
        //  res = arma::join_cols(z, res);
        //}
        // I THINK the below takes care of the corner case mentioned above....
        csy = arma::join_cols(arma::zeros<arma::vec>(1), csy);
        arma::vec res = csy.elem(bin_to_t);
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
      //arma::vec br;
      //unsigned int null_bin;
      //arma::uvec bin_to_t;
      //arma::umat t_to_bin; // either just a std::vector, or a matrix with 2 columns (first is x-values of non-null bin, second is index of non-null bin values)
      //arma::uvec order_t_low = {0};
      //arma::uvec order_t_high = {0};
  };

  // function to see if number of unique elements in a column is 1, 2, or > (return 3)
  inline unsigned int count_unique(const arma::vec& x) {
      arma::vec unique_vals(2, arma::fill::value(x[0]));
      for (size_t i = 1; i < x.size(); i++) {
          if (x[i] != unique_vals[0] && x[i] != unique_vals[1]) {
              if (unique_vals[0] != unique_vals[1]) {
                  return 3;
              }
              else {
                  unique_vals[0] = x[i];
              }
          }
      }
      if (unique_vals[0] != unique_vals[1]) {
          return 2;
      }
      return 1;
  }

  // for when X  is sparse matrix
  inline unsigned int count_unique(const arma::sp_mat& x) {
      arma::vec unique_vals(2, arma::fill::value((double)x.at(0, 0)));
      for (size_t i = 1; i < x.n_rows; i++) {
          if (x.at(i, 0) != unique_vals[0] && x.at(i, 0) != unique_vals[1]) {
              if (unique_vals[0] != unique_vals[1]) {
                  return 3;
              }
              else {
                  unique_vals[0] = x.at(i, 0);
              }
          }
      }
      if (unique_vals[0] != unique_vals[1]) {
          return 2;
      }
      return 1;
  }
  
  class StumpsMatrix {
    public:
      // constructor when using dense X
      StumpsMatrix(arma::mat& X, Rcpp::Nullable<Rcpp::IntegerVector>& _include_linear, Rcpp::Nullable<Rcpp::IntegerVector>& _include_stumps, unsigned int num_cuts, bool use_quants, unsigned int scale_X, unsigned int _ncores) {
        ncores = _ncores;
        #if defined(_OPENMP)
            omp_set_num_threads(ncores);
        #endif
        include_linear = arma::uvec(X.n_cols);
        include_stumps = arma::uvec(X.n_cols);
        bool lin_is_null = _include_linear.isNull();
        bool stumps_is_null = _include_stumps.isNull();
        if (!lin_is_null) {
            include_linear = Rcpp::as<arma::uvec>(_include_linear);
        }
        if (!stumps_is_null) {
            include_stumps = Rcpp::as<arma::uvec>(_include_stumps);
        }
        if (lin_is_null || stumps_is_null) {
            #if defined(_OPENMP)
                #pragma omp parallel for schedule(dynamic)
            #endif
            for (size_t i = 0; i < X.n_cols; i++) {
                unsigned int n_unique = count_unique(X.col(i));
                if (lin_is_null && n_unique > 1) {
                    include_linear[i] = 1;
                }
                if (stumps_is_null && n_unique > 2) {
                    include_stumps[i] = 1;
                }
            }
        }
        ncol_lin = arma::sum(include_linear); 
        arma::uvec which_incl_stumps = arma::find(include_stumps);
        //blocks = std::vector<MatrixBlock> (which_incl_stumps.size(), MatrixBlock());
        bool any_incl_lin = (ncol_lin > 0);
        if (any_incl_lin) {
          blocks = std::vector< std::shared_ptr< MatrixBlock > > (which_incl_stumps.size() + 1);
        } else {
          blocks = std::vector< std::shared_ptr< MatrixBlock > > (which_incl_stumps.size());
        }
        
        // first, make stumps blocks (in parallel)
        #if defined(_OPENMP)
          #pragma omp parallel for schedule(dynamic)
        #endif
        for (size_t i = 0; i < which_incl_stumps.size(); i++) {
          // blocks[i] = TFGMatrixBlock(X.col(which_incl_stumps[i]), cuts[which_incl_stumps[i]]);
          blocks[i + any_incl_lin].reset( new TFGMatrixBlock(X.col(which_incl_stumps[i]), num_cuts, use_quants) );
        }
        
        // now, add linear part
        // scale_X = 1, or 2 for "sd", or "max" scaling (respectively) (anything else is no scaling)
        if (any_incl_lin) {
          arma::vec col_scale_factors;
          if (scale_X == 1) {
              col_scale_factors = arma::trans(arma::stddev(X.cols(arma::find(include_linear))));
          } else if (scale_X == 2) {
              col_scale_factors = arma::trans(arma::max(arma::abs(X.cols(arma::find(include_linear)))));
          } else {
              col_scale_factors = arma::ones(ncol_lin);
          }
          blocks.front().reset( new DenseNumericMatrixBlock(X, include_linear, col_scale_factors) );
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
      StumpsMatrix(arma::sp_mat& X, Rcpp::Nullable<Rcpp::IntegerVector>& _include_linear, Rcpp::Nullable<Rcpp::IntegerVector>& _include_stumps, unsigned int num_cuts, bool use_quants, unsigned int scale_X, unsigned int _ncores) {
          ncores = _ncores;
          #if defined(_OPENMP)
            omp_set_num_threads(ncores);
          #endif
          include_linear = arma::uvec(X.n_cols);
          include_stumps = arma::uvec(X.n_cols);
          bool lin_is_null = _include_linear.isNull();
          bool stumps_is_null = _include_stumps.isNull();
          if (!lin_is_null) {
              include_linear = Rcpp::as<arma::uvec>(_include_linear);
          }
          if (!stumps_is_null) {
              include_stumps = Rcpp::as<arma::uvec>(_include_stumps);
          } 
          if (lin_is_null || stumps_is_null) {
              #if defined(_OPENMP)
                #pragma omp parallel for schedule(dynamic)
              #endif
              for (size_t i = 0; i < X.n_cols; i++) {
                  unsigned int n_unique = count_unique(X.col(i));
                  if (lin_is_null && n_unique > 1) {
                      include_linear[i] = 1;
                  }
                  if (stumps_is_null && n_unique > 2) {
                      include_stumps[i] = 1;
                  }
              }
          }
          ncol_lin = arma::sum(include_linear);
          arma::uvec which_incl_stumps = arma::find(include_stumps);
        //blocks = std::vector<MatrixBlock> (which_incl_stumps.size(), MatrixBlock());
        bool any_incl_lin = (ncol_lin > 0);
        if (any_incl_lin) {
          blocks = std::vector< std::shared_ptr< MatrixBlock > > (which_incl_stumps.size() + 1);
        } else {
          blocks = std::vector< std::shared_ptr< MatrixBlock > > (which_incl_stumps.size());
        }
        
        // first, make stumps blocks (in parallel)
        #if defined(_OPENMP)
          #pragma omp parallel for schedule(dynamic)
        #endif
        for (size_t i = 0; i < which_incl_stumps.size(); i++) {
          // blocks[i] = TFGMatrixBlock(X.col(which_incl_stumps[i]), cuts[which_incl_stumps[i]]);
          blocks[i + any_incl_lin].reset( new TFGMatrixBlock(X.col(which_incl_stumps[i]), num_cuts, use_quants) );
        }
        
        // now, add linear part
        // scale_X = 1, or 2 for "sd", or "max" scaling (respectively) (anything else is no scaling)
        if (any_incl_lin) {
            arma::vec col_scale_factors;
            if (scale_X == 1) {
                col_scale_factors = arma::trans(arma::stddev(arma::mat(X.cols(arma::find(include_linear)))));
            } else if (scale_X == 2) {
                col_scale_factors = arma::trans(arma::max(arma::abs(X.cols(arma::find(include_linear)))));
            } else {
                col_scale_factors = arma::ones(ncol_lin);
            }
            blocks.front().reset(new SparseNumericMatrixBlock(X, include_linear, col_scale_factors));
        }
        
        nrow = X.n_rows;
        ncol = 0;
        // get ncol, and col_offsets for each block
        for (size_t i = 0; i < blocks.size(); i++) {
          blocks[i].get()->col_offset = ncol;
          ncol += blocks[i].get()->ncol;
        }
      }

      // constructor when input is StumpsMatrix made with dense X
      StumpsMatrix(arma::mat& X, Rcpp::XPtr<StumpsMatrix> SM_in) {
          ncores = SM_in.get()->ncores;
          #if defined(_OPENMP)
            omp_set_num_threads(ncores);
          #endif   
          include_linear = SM_in.get()->include_linear;
          include_stumps = SM_in.get()->include_stumps;
          ncol_lin = SM_in.get()->ncol_lin;
          arma::uvec which_incl_stumps = arma::find(include_stumps);
          //blocks = std::vector<MatrixBlock> (which_incl_stumps.size(), MatrixBlock());
          bool any_incl_lin = (ncol_lin > 0);
          if (any_incl_lin) {
              blocks = std::vector< std::shared_ptr< MatrixBlock > >(which_incl_stumps.size() + 1);
          }
          else {
              blocks = std::vector< std::shared_ptr< MatrixBlock > >(which_incl_stumps.size());
          }

          // first, make stumps blocks (in parallel)
          #if defined(_OPENMP)
            #pragma omp parallel for schedule(dynamic)
          #endif
          for (size_t i = 0; i < which_incl_stumps.size(); i++) {
              // blocks[i] = TFGMatrixBlock(X.col(which_incl_stumps[i]), cuts[which_incl_stumps[i]]);
              blocks[i + any_incl_lin].reset(new TFGMatrixBlock(X.col(which_incl_stumps[i]), SM_in.get()->blocks[i + any_incl_lin]));
          }

          // now, add linear part
          // scale_X = 1, or 2 for "sd", or "max" scaling (respectively) (anything else is no scaling)
          if (any_incl_lin) {
              blocks.front().reset(new DenseNumericMatrixBlock(X, SM_in.get()->blocks[0]));
          }

          nrow = X.n_rows;
          ncol = 0;
          // get ncol, and col_offsets for each block
          for (size_t i = 0; i < blocks.size(); i++) {
              blocks[i].get()->col_offset = ncol;
              ncol += blocks[i].get()->ncol;
          }
      }

      // constructor when input is StumpsMatrix made with sparse X
      StumpsMatrix(arma::sp_mat& X, Rcpp::XPtr<StumpsMatrix> SM_in) {
          ncores = SM_in.get()->ncores;
          #if defined(_OPENMP)
            omp_set_num_threads(ncores);
          #endif   
          include_linear = SM_in.get()->include_linear;
          include_stumps = SM_in.get()->include_stumps;
          ncol_lin = SM_in.get()->ncol_lin;
          arma::uvec which_incl_stumps = arma::find(include_stumps);
          //blocks = std::vector<MatrixBlock> (which_incl_stumps.size(), MatrixBlock());
          bool any_incl_lin = (ncol_lin > 0);
          if (any_incl_lin) {
              blocks = std::vector< std::shared_ptr< MatrixBlock > >(which_incl_stumps.size() + 1);
          }
          else {
              blocks = std::vector< std::shared_ptr< MatrixBlock > >(which_incl_stumps.size());
          }

          // first, make stumps blocks (in parallel)
          #if defined(_OPENMP)
            #pragma omp parallel for schedule(dynamic)
          #endif
          for (size_t i = 0; i < which_incl_stumps.size(); i++) {
              // blocks[i] = TFGMatrixBlock(X.col(which_incl_stumps[i]), cuts[which_incl_stumps[i]]);
              blocks[i + any_incl_lin].reset(new TFGMatrixBlock(X.col(which_incl_stumps[i]), SM_in.get()->blocks[i + any_incl_lin]));
          }

          // now, add linear part
          // scale_X = 1, or 2 for "sd", or "max" scaling (respectively) (anything else is no scaling)
          if (any_incl_lin) {
              blocks.front().reset(new SparseNumericMatrixBlock(X, SM_in.get()->blocks[0]));
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
          #pragma omp parallel for reduction(+:Xb) schedule(dynamic)
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
          #pragma omp parallel for schedule(dynamic)
        #endif
        for(size_t i = 0; i < blocks.size(); i++) {
          Xty.rows(blocks[i].get()->col_offset, blocks[i].get()->col_offset + blocks[i].get()->ncol - 1) = blocks[i].get()->compute_Xty(y, X_avg.rows(blocks[i].get()->col_offset, blocks[i].get()->col_offset + blocks[i].get()->ncol - 1));
        }
        return(Xty);
      }
      
      arma::vec compute_X2b(const arma::vec& b, const arma::vec& X_avg) {
        arma::vec X2b(nrow, arma::fill::zeros);
        #if defined(_OPENMP)
          #pragma omp parallel for reduction(+:X2b) schedule(dynamic)
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
          #pragma omp parallel for schedule(dynamic)
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
      unsigned int ncol_lin;
      
    protected:
      std::vector< std::shared_ptr< MatrixBlock > > blocks;
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

