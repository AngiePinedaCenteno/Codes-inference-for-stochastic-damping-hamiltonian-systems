# Codes-inference-for-stochastic-damping-hamiltonian-systems
This repository contain the codes used for the simulations
of our article: 

An overlook on statistical inference issues for
stochastic damping hamiltonian systems under
the fluctuation-dissipation condition
with some modifications that makes the
codes show  the results much  more quickly.

In the folder named CodesFluctuationDissipation you will
find the following functions of matlab:

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
1)  Execute the function named Inputs_FluctuationDissipation.m
     This function generates the input variables of the function named                   
     Outputs_FluctuationDissipation.m which is the principal function.

     This function  (Inputs_FluctuationDissipation.m ), don't has  input
     variables. If you want to modify the variables of the simulation you
     must modify the values in the script of the function
    Inputs_FluctuationDissipation.m .
     Below is a view of this  function. At right of the script are the descrip-
     tion of the variables.

2)  Execute the function named Outputs_FluctuationDissipation.m
     This function generates all the estimations corresponding to the system
     of  SDE that we are studing.

The rest of the functions in the folder are  necessary to the execution of the principal
function named Outputs_FluctuationDissipation.m. All the functions contain its
description at the top of its script.

Finally, to obtain  the comparative graphics of the theoretical functions versus its
estimates, you must to execute the function named FigModule2pi.m. To do this you
only must to  replace the input variable casG for the number 1, 2, 3, 4 or 5.  For example
if you make casG=1, you will obtain the   graphics of the densities and the bidimensional
graphics  of the drift an the derivative. The full description is at the top of the script of the
 function.
