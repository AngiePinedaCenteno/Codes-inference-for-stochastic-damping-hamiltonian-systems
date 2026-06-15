function [UX,UY,p,px,py,dKp,dKpx,dKpy,Dp,DdKp,g,driftE,driftE0And1,Hn,DD,MbetgV,EbetgVM,xd,yd,sxd,syd,hn,ElapsedTime]=Outputs_FluctuationDissipation(n,gamm,rm,bet,x0,y0,dx,dy,alpha1,alpha2,varep,period,casV,casS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes non parametric estimates of the invariant density 
% of the process (X,Y)=(X_t,Y_t)_(t>=0) which resolves the system:   
%                  dX_t=Y_tdt  
%                  dY_t=sqrt(2/bet)*s*dW_t -(s^2*Y_t+gV(X_t))*dt,
% where: bet is a constant, W is a Brownian motion and gV denotes the 
% gradient of the potencial V.
% Also,this function computes non parametric estimations of the drift term, 
% the potential and the variance of the system defined above.
%
% The input variables of this function are described at the top of the 
% function: Inputs_FluctuationDissipation.m
% The outputs of this function are:
% .- UX and UY denote the approximated discrete subsampling of the trajec-
%    tories for the position and the velocity that resolve the above SDE
% .- p is a discretization of the invariant density
% .- px and py are the corresponding marginal invariant densities in x and
%    y direction
% .- Dp denote the discretization of the derivative with respect to x of
%    the invariant density. 
% .- g denote a discretization of the drift
% .- DD is a discretization of the drift multiplied by the invariant
%    density.
% .- MbetgV is a discretization of the gradient of the potential multiplied
%    by the constant -bet.
% .- dKp denotes the estimation of the invariant density
% .- dKpx and dKpy are the estimates of the marginal invariant densities in
%    x and y direction
% .- DdKp denotes the estimation of the derivative with respect to x of the 
%    invariant density. 
% .- driftE denotes the estimate of the drift
% .- driftE0And1 stores in the first row the drift estimate for y=0, and in
%    the second row stores the drift estimate for y=1.
% .- Hn is the estimate of the drift multiplied by the invariant density.
% .- EbetgV is the estimation of the gradient of the potential multiplied
%    by the constant -bet.
% .- xd is the discretization of the domain of the dendity in x-direction
% .- yd is the discretization of the domain of the dendity in y-direction
% .- sxd and syd are the lengths of the vectors xd and yd respectively
% .- ElapsedTim is the time in minutes of the execution  of this  function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic                                                                        %% starts a stopwatch timer to measure performance of this fnction (Outputs_FluctuationDissipation.m)
%% Parameters
b1n=n^(-alpha1); b2n=n^(-alpha2);                                          %% b1n and b2n are the bandwidths of the non parametric estimator of the density p_s^{\bet}                                                       
%% Generation of te trajectories of the position and velocity of the system
[UX,UY,hn]=SDE(n,gamm,rm,bet,x0,y0,casV,casS);                             %% This function generates the approximation of a trajectory of X an a trajector for Y. 
if ((casV==1) || (casV==2) )
    %% Transformation modulo 2*pi centered
    UX=UX+(period/2);
    UX=UX-(floor((1/period)*UX)*period);
    UX=UX-(period/2);
end
%% DensityAndDrift.m generates the estimate of: the invariant density, the 
% derivative wih respect to x of the invariant density,and the drisft.
[dKp,dKpx,dKpy,driftE,DdKp,Hn,xd,yd,sxd,syd]=DensityDerivativeDrift(UX,UY,b1n,b2n,hn,dx,dy,n,varep);
%% DensityAndDrift_0And1.m generates the estimate of the drisft of the SDE 
% in y=0 and y=1.
[driftE0And1]=DensityAndDrift_0And1(UX,UY,xd,sxd,n,b1n,b2n,hn,varep) ;
%% Theoretical densities
[p,px,py,Dp,g,DD,MbetgV] = TheoreticalFunctions(xd,yd,bet,dx,dy,casV,casS);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EbetgVM=DdKp./(dKp+varep);
toc                                                                        %% returns the elapsed time in seconds of this function (Outputs_FluctuationDissipation.m)
ElapsedTime=toc/60;                                                        %% time in minutes of the exection  of this  function (Outputs_FluctuationDissipation.m)
end