# 1. loglike function in C++  
The cmake file right now test the origianl heston model marginal log-likelihood. (setting jump parameters to 0) It uses the old data 'YObs' from matlab.

## Using Cmake, cross-platform
First, check if c++ boost library is installed.  

```{bash}
rm -rf build
mkdir build && cd build
cmake ..
make
./a.out
```

## generate shared/static library without using CMake, on Linux or Mac OS X
#### shared
```{bash}
g++ --std=c++11 test_loglike.cpp -Iinclude -fPIC -shared -o libtestloglike.so
g++ --std=c++11 test_loglike.cpp -Iinclude -L. -ltestloglike -o a.out
./a.out
```
#### static
```{bash}
g++ -c --std=c++11 test_loglike.cpp -Iinclude
ar -crv libloglike.a test_loglike.o
g++ --std=c++11 test_loglike.cpp -Iinclude -L. -lloglike -o a.out 
./a.out
```

## For VS Studio in Windows using CMake  
Simply open the folder containing cmake.  
Details:`https://blogs.msdn.microsoft.com/vcblog/2016/10/05/cmake-support-in-visual-studio/`

# 2. Optimization in R using package 'nloptr'  
See the 'R' folder.

