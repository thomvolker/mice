#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

//' Find index of matched donor units
//'
//' @param yhatobs   Numeric vector with values from donor cases.
//' @param yhatmis   Numeric vector with values from target cases.
//' @param k   Integer, number of unique donors from which a random draw is made.
//' For \code{k = 1} the function returns the index in \code{d} corresponding
//' to the closest unit. For multiple imputation, the
//' advice is to set values in the range of \code{k = 5} to \code{k = 10}.
//' @param cond    Vector with values of the same length as yhatmis, with a
//' value to exclude from the donor pool for each missing observation.
//' @param ytrue   The observed y-values from which the donors are drawn, to
//' discard the observations from this set.
//' @return An integer vector with \code{length(t)} elements. Each
//' element is an index in the array \code{d}.
//' @details
//' For each element in \code{t}, the method finds the \code{k} nearest
//' neighbours in \code{d}, randomly draws one of these neighbours, and
//' returns its position in vector \code{d}.
//'
//' Fast predictive mean matching algorithm in seven steps:
//'
//' 1. Shuffle records to remove effects of ties
//'
//' 2. Obtain sorting order on shuffled data
//'
//' 3. Calculate index on input data and sort it
//'
//' 4. Pre-sample vector \code{h} with values between 1 and \code{k}
//'
//' For each of the \code{nmis} elements in \code{t}:
//'
//'   5. find the two adjacent neighbours
//'
//'   6. find the \code{h_i}'th nearest neighbour
//'
//'   7. store the index of that neighbour
//'
//' Return vector of \code{nmis} positions in \code{d}.
//'
//' We may use the function to perform predictive mean matching under a given
//' predictive model. To do so, specify both \code{d} and \code{t} as
//' predictions from the same model. Suppose that \code{y} contains the observed
//' outcomes of the donor cases (in the same sequence as \code{d}), then
//' \code{y[matchindex(d, t)]} returns one matched outcome for every
//' target case.
//'
//' See \url{https://github.com/amices/mice/issues/236}.
//' This function is a replacement for the \code{matcher()} function that has
//' been in default in \code{mice} since version \code{2.22} (June 2014).
//' @examples
//' set.seed(1)
//'
//' # Inputs need not be sorted
//' d <- c(-5, 5, 0, 10, 12)
//' t <- c(-6, -4, 0, 2, 4, -2, 6)
//'
//' # Index (in vector a) of closest match
//' idx <- matchindex(d, t, 1)
//' idx
//'
//' # To check: show values of closest match
//'
//' # Random draw among indices of the 5 closest predictors
//' matchindex(d, t)
//'
//' # An example
//' train <- mtcars[1:20, ]
//' test <- mtcars[21:32, ]
//' fit <- lm(mpg ~ disp + cyl, data = train)
//' d <- fitted.values(fit)
//' t <- predict(fit, newdata = test)  # note: not using mpg
//' idx <- matchindex(d, t)
//'
//' # Borrow values from train to produce 12 synthetic values for mpg in test.
//' # Synthetic values are plausible values that could have been observed if
//' # they had been measured.
//' train$mpg[idx]
//'
//' # Exercise: Create a distribution of 1000 plausible values for each of the
//' # twelve mpg entries in test, and count how many times the true value
//' # (which we know here) is located within the inter-quartile range of each
//' # distribution. Is your count anywhere close to 500? Why? Why not?
//' @author Stef van Buuren, Nasinski Maciej, Alexander Robitzsch

//' @export 
// [[Rcpp::export]]
arma::vec matchindex_dev(arma::vec yhatobs, arma::vec yhatmis, int k, arma::vec cond, arma::vec ytrue) {
  
  Environment base("package:base");
  Function sample = base["sample"];
  
  // declarations
  int nobs = yhatobs.size();
  int nmis = yhatmis.size();
  
  arma::uvec ishuf = as<arma::uvec>(sample(nobs));
  ishuf = ishuf - 1;
  arma::vec yshuf = yhatobs(ishuf);
  arma::uvec id = arma::sort_index(yshuf);
  arma::vec ysort = yshuf(id);
  
  // 4. Pre-sample nmis values between 1 and k
  // restrict 1 <= k <= nobs
  if (k < 1) {k = 1;}
  if (k > nobs) {k = nobs;}
  // k = arma::min(arma::uvec({max(arma::uvec({1, k})), nobs}));
  arma::uvec h = as<arma::uvec>(sample(k, nmis, Rcpp::_["replace"] = true));
  
  arma::vec idx = arma::zeros<arma::vec>(nmis);
  
  for (int i = 0; i < nmis; i++) {
    double val = yhatmis(i);
    int ex = cond(i);
    arma::vec yuse = ysort(arma::find(ytrue(id) != ex));
    arma::uvec iduse = id(arma::find(ytrue(id) != ex));
    int nuse = yuse.size();
    int hi = h(i);
    int count = 0;
    
    // 5. find the two adjacent neighbours
    arma::vec::iterator iter = std::lower_bound(yuse.begin(), yuse.end(), val);
    int r = iter - yuse.begin();
    int l = r - 1;
    
    while (count < hi && l >= 0 && r < nuse) {
      if (val - yuse(l) < yuse(r) - val) {
        idx(i) = iduse(l--);
      } else {
        idx(i) = iduse(r++);
      }
      count++;
    }
    
    // If right side is exhausted, take left elements
    while (count < hi && l >= 0) {
      idx(i) = iduse(l--);
      count++;
    }
    
    // If left side is exhausted, take right elements
    while (count < hi && r < nuse) {
      idx(i) = iduse(r++);
      count++;
    }
  }
  return idx + 1;
}
 