# C++ part  
I used CMake for the project. `include` contain the header-only implementation of the functions. `test` contains some cpp files to test these functions. Specially,
the cmake file right now test `test_loglike.cpp`, which check the loglikelihood results for the original Heston model.(setting jump parameters to 0) It uses the Jiayun's parameters as well as path of 'YObs'.

### Compile instruction 
#### Using Cmake
First, check if c++ boost library is installed. (used for the incomplete gamma function).  
Then from commandline:  
```{bash}
rm -rf build
mkdir build && cd build
cmake ..
make
./a.out
```
Another way might be:  
### generate shared/static library
##### shared
```{bash}
g++ --std=c++11 test_loglike.cpp -Iinclude -fPIC -shared -o libtestloglike.so
g++ --std=c++11 test_loglike.cpp -Iinclude -L. -ltestloglike -o a.out
./a.out
```
##### static
```{bash}
g++ -c --std=c++11 test_loglike.cpp -Iinclude
ar -crv libloglike.a test_loglike.o
g++ --std=c++11 test_loglike.cpp -Iinclude -L. -lloglike -o a.out 
./a.out
```

### For VS Studio in Windows using CMake  
Simply open the folder containing cmake.  
Details:`https://blogs.msdn.microsoft.com/vcblog/2016/10/05/cmake-support-in-visual-studio/`

# R
The R script need some libraries. Install them by typing in the R console:

```
install.packages(c("Rcpp", "BH", "ESGtoolkit", "nloptr"))
```

`loglike_port_r.cpp` port the C++ function into R.    
Other comments are made within `heston_svj_mle_complete.Rmd` file.  

