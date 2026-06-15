function [ n,gamm,rm,bet,x0,y0,dx,dy,alpha1,alpha2,varep,period,casV,casS ] = Inputs_FluctuationDissipation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function determine the inputs variables of the simulations, i.e. the
% parameters of the SDE:
%                  dX_t=Y_tdt  
%                  dY_t=sqrt(2/bet)*s*dW_t -(s^2*Y_t+gV(X_t))*dt,
% and the variables of the density estimator 'dKp',the derivative with res-
% pect to x of the estimator of the invariant density 'Ddkp', a function of 
% the potential estimator ''EbetgVM' and the drift estimator 'driftE'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=10^4;                                                                    %% n is the number of time steps for the approximations of the trajectories X and Y of the solution of the SDE.
gamm=0.15;                                                                 %% gamm is the potency of n to form the value of h_n (h_n=n^(-gamm))
rm=30;                                                                     %% rm determines the relation delta=hn/rm between the time steps delta and hn 
bet=2;                                                                     %% bet is the covnstant bet in the difussion term of the above system
x0=0;                                                                      %% x0 is the inital condition of the system in the x direction
y0=0;                                                                      %% y0 is the inital condition of the system in the y direction
alpha1=0.2;                                                                %% alpha1 is the potency of n to form the value of b_1n (b_1n=n^(-alpha1)) 
alpha2=0.2;                                                                %% alph2 is the potency of n to form the value of b_2n (b_2n=n^(-alpha2))                                     
dx=1/15;                                                                   %% dx is the grid spacing in the x direction used in the discretization of the domain of the invariant density                                                      
dy=1/15;                                                                   %% dy is the grid spacing in the y direction used in the discretization of the domain of the invariant density
varep=0.0001;                                                              %% varep is a little valued that is added to the estimate of the invariant density to make it nonzero.                    
period=2*pi;                                                               %% period is  the value of the period of the potential function thar defines the SDE
casV=2;                                                                    %% casV is the variable which value determines the potential that will be used to define the system os SDEs.
casS=1;                                                                    %% casS is the variable which value determines the variance that will be used to define the system os SDEs.
end

