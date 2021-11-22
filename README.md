# Multivariate-student-t-parameter-estimate
This function returns the approximate parameters from a n-by-p sample assuming to be a p dimentional multivariate t distribution, by method Batch Approximation Algorithm attributed to Aeschliman, Park & Kak(2010)

Input

- T is a n-by-p sample matrix with n observations and p variables

- knownpara1 and knownpara2 (Optional) are the known parameters (mu,sacle or degree of freedom) which force the estimate to be the known true ones. See below for more information about Options

Output

- muhat returns the approximate p-by-1 location vector

- Sigmahat returns the approximate p-by-p scale matrix

- nuhat returns the approximate degree of freedom

Optional Inputs (knownpara1 & knownpara2)

- known mu must be p-by-1 vector

- known scale matrix must be p-by-p matrix, when you only know the scale matrix is a scaled diagonal matrix, specify the 'knownpara' Option to be zeros(p)
      
- known degree of freedom must be a 1-by-1 scaler
