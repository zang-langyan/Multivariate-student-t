# Multivariate-student-t
This repository contain functions based on Matlab to compute the Multivariate t distribution and conditional density and parameters.

MVT.m can only compute the Bivariate t distribution density for now

MVT_Con computes the Conditional distribution parameters X2|X1 based on Peng Ding(2016)[1] theory work on conditional distribution of multivariate t

MVTrand generates the random numbers which are multivariate t distributed

MLE_con_t estimates the Conditional X2|X1 parameters via Maximum likelihood Estimator and titer.m attributed to Prof. Paolella

MVTpara returns the approximate parameters from a n-by-p sample assuming to be a p dimentional multivariate t distribution, by method Batch Approximation Algorithm attributed to Aeschliman, Park & Kak(2010)[2]

# References
[1] Peng Ding. “On the Conditional Distribution of the Multivariate t Distribution”. In: The American Statistician 70.3 (July 2, 2016), pp. 293–295. issn: 0003-1305, 1537-2731. doi: 10.1080/00031305.2016.1164756. url: https://www.tandfonline.com/doi/full/ 10.1080/00031305.2016.1164756 (visited on 11/22/2021).

[2] David Hutchison et al. “A Novel Parameter Estimation Algorithm for the Multivariate t-Distribution and Its Application to Computer Vision”. In: Computer Vision – ECCV 2010. Ed. by Kostas Daniilidis, Petros Maragos, and Nikos Paragios. Vol. 6312. Series Title: Lecture Notes in Computer Science. Berlin, Heidelberg: Springer Berlin Heidelberg, 2010, pp. 594–607. isbn: 978-3-642-15551-2 978-3-642-15552-9. doi: 10.1007/978-3-642- 15552-9_43. url: http://link.springer.com/10.1007/978-3-642-15552-9_43 (visited on 11/22/2021).
