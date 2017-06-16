# ConcurrentLab_3_4

This is a CodeBlocks project, therefore the easiest way to access the code will be to open `ConcurrentProg_Lab3_4.cbp`(CodeBlocks project file) in CodeBlocks.

Alternatively, you can run `main.cpp` 

# Matrix class

<b>Initializing:</b> 

`Matrix a(5);` : initializes a square matrix of size 5


<br><b>Size</b>

`a.getSize();` : returns 5


<br><b>Fill with random values</b>

`a.generateRandomValues();` : fills the matrix with random values


<br><b>Display elements:</b> 

`a.displayValues();` display all elements in the matrix


<br><b>Accessing elements:</b> 

`a(0,0);` : first element in matrix

`a(4,4);` : last element in matrix 

`a(1,4);` : 2nd row 5th element



<br>

# Functions


<br><b>Sequential multiplication:</b> 

`Matrix C = seq_mat_mul(Matrix A, Matrix B)` : calculates sequential matrix multiplication of A and B and returns a new `Matrix`


<br><b>Parallel multiplication:</b> 

`Matrix C = parallel_mat_mul(Matrix A, Matrix B)` : calculates parallel for loop matrix multiplication of A and B using OpenMP and returns a new `Matrix`
