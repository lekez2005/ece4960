### Part I: Ground Truth 
Verify original full matrix prints correctly as a sparse matrix
|1 | 2 | 0 | 0 | 3 |
|4 | 5 | 6 | 0 | 0 |
|0 | 7 | 8 | 0 | 9 |
|0 | 0 | 0 | 10 | 0 |
|11 | 0 | 0 | 0 | 12 |
Verify norm works correctly using slightly modified test matrix
|1 | 2 | 2.5 | 0 | 3 |
|4 | 5 | 6 | 0 | 0 |
|0 | 7 | 8 | 0 | 9 |
|0 | 0 | 0 | 10 | 0 |
|11 | 0 | 0 | 0 | 12 |
Expected norm is 6.25, Actual norm is 6.25 

Result of permute (0, 2): 
|0 | 7 | 8 | 0 | 9 |
|4 | 5 | 6 | 0 | 0 |
|1 | 2 | 0 | 0 | 3 |
|0 | 0 | 0 | 10 | 0 |
|11 | 0 | 0 | 0 | 12 |
Norm of sparse and full matrix solutions: 	 0

Result of permute (0, 4): 
|11 | 0 | 0 | 0 | 12 |
|4 | 5 | 6 | 0 | 0 |
|1 | 2 | 0 | 0 | 3 |
|0 | 0 | 0 | 10 | 0 |
|0 | 7 | 8 | 0 | 9 |
Norm of sparse and full matrix solutions: 	 0

Result of 3.0*row[0] + row[3] 
|1 | 2 | 0 | 0 | 3 |
|4 | 5 | 6 | 0 | 0 |
|0 | 7 | 8 | 0 | 9 |
|3 | 6 | 0 | 10 | 9 |
|11 | 0 | 0 | 0 | 12 |
Norm of sparse and full matrix solutions: 	 0

Result of -4.4*row[4] + row[1] 
|1 | 2 | 0 | 0 | 3 |
|-44.4 | 5 | 6 | 0 | -52.8 |
|0 | 7 | 8 | 0 | 9 |
|0 | 0 | 0 | 10 | 0 |
|11 | 0 | 0 | 0 | 12 |
Norm of sparse and full matrix solutions: 	 0

For x = 
5, 	4, 	3, 	2, 	1
Result of A*x = 
16, 	58, 	61, 	20, 	67
Norm of sparse and full matrix solutions: 	 0
### Part 2 
Sum using product = 101.594, Sum using values = 101.594
Residue = 8.52651e-13
Time elapsed = 364.198 ms 
Peak Memory usage using /proc/self/status -> VmPeak :
VmPeak:	   18220 kB
