# SMINT
The clustering technique to cluster non-Gaussian functional data without the prior information of the number of clusters.

simulation 1.R—simulation 5.R are the the main code corresponding to the simulations in the paper, real.R is the main code for ADNI data. It uses functions defined in the files EM_estimate1.R, EM_estimate2.R, EM_estimate4.R, original function.R, first_der function.R, H_est.R, identify.R.

EM_estimate: the function for estimating Omega_n given H

H_est: the function for estimating H  given Omega_n

# The arguments in the main code file:

The number of datasets: N

sample size: n

penalty parameter: lambda

the mixing probabilities: PI

the number of eigenfunctions: K_g

the true number of clusters: true_num

the varainces of errors: ggamma

the transformation function: H 

the inverse function of H: SQ

the mean functions:  g1f, g2f, g3f

the eigenfunctions: eigenfun

the initial number of clusters: num

the true clutser id: cluster_result[,1:n]

the functional observation：YY, a list

the obvserved time: TT, a list




# The estimates in the main code file:


the mixing probabilities: prob, N * true_num matrix

the basis coefficient of mean functions: parmet, N * q_n * true_num  array

the varainces of errors: stderror, N * true_num  matrix

the eigenvalue: eigenvalue, N * K_g * K_g * true_num  array

the basis coefficient of eigenfunctions: eigencoe, N * K_g * q_n * true_num  array

the value of H: H_fun, N * 170 matrix

the number of clusters: num_group, N-dimensional vector

the estimated clutser id: cluster_result[,(n+1):(2 * n)], N * n matrix

