# Codes-inference-for-stochastic-damping-hamiltonian-systems
This repository contain the codes I made for the numerical experiments of one of my researh collaborations, the article :

P. Cattiaux, J. R. León, A. A. Pineda Centeno & C. Prieur (2017). An overlook
on statistical inference issues for stochastic damping hamiltonian systems under the fluctuationdissipation
condition, Statistics, 51:1, 11-29, DOI: 10.1080/02331888.2016.1259807

In this repository you will find the following functions that I created by using matlab :

.- Inputs_FluctuationDissipation.m

.- Outputs_FluctuationDissipation.m

.- SDE.m

.- DensityDerivativeDrift.m

.- DensityAndDrift_0And1.m

.- TheoreticalFunctions.m

.- EK.m

.- Pot.m

.- S.m

To obtain the data set of simulations you must :

1)  Execute the function named Inputs_FluctuationDissipation.m. 
     This function generates the input variables of the function named                   
     Outputs_FluctuationDissipation.m which is the principal function.

2)  Execute the function named Outputs_FluctuationDissipation.m. 
     This function generates all the estimations corresponding to the system
     of  SDE we studied in the article.

The rest of the functions in the folder are  necessary to the execution of the principal
function named Outputs_FluctuationDissipation.m. All the functions contain its
description at the top of its script.

Finally, to obtain  the comparative graphics of the theoretical functions versus its
estimates, you must to execute the function named FigModule2pi.m. 
