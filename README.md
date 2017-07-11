# ConcurrentLab_3_4

The latest working code is in <b>mat_mul_v4</b> and the compiler used is GNU GCC (https://gcc.gnu.org/)

Note: OpenMP must be available in the compiler. (OpenMP Compilers: http://www.openmp.org/resources/openmp-compilers/)

# Sequential

To compile: 
<pre>g++ -std=gnu++11 -o sequential sequential.cpp</pre>

To run: 
<pre>sequantial.exe</pre>

# Parallel - not optimized

To compile:
<pre>g++ -std=gnu++11 -g -Wall -fopenmp -o parallel_not_opt parallel_not_opt.cpp</pre>

To run:
<pre>parallel_not_opt</pre>

# Parallel - optimized

To compile:
<pre>g++ -std=gnu++11 -g -Wall -fopenmp -mavx -o parallel_opt parallel_opt.cpp</pre>

To run:
<pre>parallel_opt 128</pre> <i>- where 128 is the L2 cache size</i>
