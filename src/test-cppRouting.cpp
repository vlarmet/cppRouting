#include <testthat.h>
#include <Rcpp.h>
#include "dijkstra_mat.h"

int twoPlusTwo();
//Rcpp::NumericMatrix Dijkstra_mat(std::vector<int> gfrom,std::vector<int> gto,std::vector<double> gw,int NbNodes,std::vector<int> dep, std::vector<int> arr);
// Normally this would be a function from your package's
// compiled library -- you might instead just include a header
// file providing the definition, and let R CMD INSTALL
// handle building and linking.

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.
context("Sample unit tests") {

  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  test_that("two plus two equals four") {
    expect_true(twoPlusTwo() == 4);
  }

}
// [[Rcpp::export]]
context("First cppRouting function") {
  test_that("function executes") {

    std::vector<int> gfrom;
    std::vector<int> gto;
    std::vector<double> gw;
    int NbNodes;
    std::vector<int> dep;
    std::vector<int> arr;
    gfrom = std::vector<int>{ 0, 0, 1, 1, 2, 2, 3, 4, 4 };
    gto = std::vector<int>{ 1, 3, 2, 4, 4, 5, 1, 3, 5 };
    gw = std::vector<double>{ 9, 2, 11, 3, 5, 12, 4, 1, 6 };
    NbNodes = 6;
    dep  = std::vector<int>{ 0, 1, 2, 3,  4, 5 };
    arr = std::vector<int>{ 0, 1, 2,  3, 4, 5 };
    Rcpp::NumericVector w = { 1, 2, 3, 4 };
    Rcpp::NumericVector e = { 0, R_NaReal, R_NaReal, R_NaReal, R_NaReal, R_NaReal, R_NaReal, 6,  0, 10,  4,  5, 0, 17, 11,  0, 15, 16, 0, 2,  4,  6,  0,  1, 0, 9,  3,  5,  7,  0, 0, 15,  9, 11, 13,  6,  0 };
    //Rcpp::NumericVector e = { 0, Rcpp::NA, Rcpp::NA, Rcpp::NA, Rcpp::NA, Rcpp::NA, 6,  0, 10,  4,  5, Rcpp::NA, 17, 11,  0, 15, 16, Rcpp::NA, 2,  4,  6,  0,  1, Rcpp::NA, 9,  3,  5,  7,  0, Rcpp::NA, 15,  9, 11, 13,  6,  0 };
    e.attr("dim") = Rcpp::Dimension(6, 6);
helloI've     Rcpp::NumericMatrix expected = Rcpp::as<Rcpp::NumericMatrix>(e);
    Rcpp::NumericMatrix output( 6, 6 );
    output = InternalDijkstra_mat(gfrom, gto, gw, NbNodes, dep, arr);
    //INFO("NbNodes" << NbNodes << "\n");
    expect_true(output.nrow() == 6);
    for (int i = 0; i < 6; i++) {
      for (int j = 0; j < 6; j++) {
        expect_true(output(i, j) == expected(i, j));
          }}
  }
}
