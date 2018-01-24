## Simple MC Pricer for Heston Model

#### Output
```
> mkdir build
> cd build
> cmake ..
> make
> ./a.out
Calculating path 10000 of 100000
Calculating path 20000 of 100000
Calculating path 30000 of 100000
Calculating path 40000 of 100000
Calculating path 50000 of 100000
Calculating path 60000 of 100000
Calculating path 70000 of 100000
Calculating path 80000 of 100000
Calculating path 90000 of 100000
Calculating path 100000 of 100000
Option Price: 6.80981
```

#### Benchmark
The exact option price is 6.8061, as reported by *Broadie and Kaya*, Exact Simulation of Stochastic Volatility and other Affine Jump Diffusion Processes
