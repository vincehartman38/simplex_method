# simplex_method
Python implementation of the two phase simplex method within linear programming using Bland's rule

The algorithm will provide one of three solutions - a bounded solution, no solution, or an infinite (unbounded) solution.

The program takes the following variables:  
int n: number of inequality equations  
int m: total variables  
list a: coefficients of inequality equations in a list of list  
list b: maximums for each of the inqualities  
list c: coefficients of the optimization function  

Example Problem:  
x + y - 3z <= 10  
5x + 10y <= -50  
3x - 2y -4z <= 9  
Maximize: -x - 6y - 3z  

You would enter the following within the input prompt:  
3 3  
1 1 -3  
-5 10 0  
3 -2 -4  
10 -50 9  
-1 -6 -3  

Bounded Solution: 10.000 0.000 5.250
