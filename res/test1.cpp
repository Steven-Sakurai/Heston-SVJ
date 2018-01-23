#include<Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
void test_copy(){
    NumericVector A = NumericVector::create(1, 2, 3);
    NumericVector B = A;
    
    Rcout << "Before: " << std::endl << "A: " << A << std::endl << "B: " << B << std::endl; 
    
    A[1] = 5; // 2 -> 5
    
    Rcout << "After: " << std::endl << "A: " << A << std::endl << "B: " << B << std::endl; 
}