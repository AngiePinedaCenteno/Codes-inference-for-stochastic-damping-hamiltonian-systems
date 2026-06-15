function [p,px,py,Dp,g,DD,MbetgV] = TheoreticalFunctions(xd,yd,bet,dx,dy,casV,casS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function compute the discretization of the invariant density, the
% drift, a function of the potential, and the variance.
% Input variables:
% .- xd is a discretization of the domain of the dendity in the x-direction
% .- yd is a discretization of the domain of the dendity in the y-direction
% The rest of the input of this function are described in the function
% Inputs_FluctuationDissipation.m 
% The outputs of this function are decribes in the function
% Outputs_FluctuationDissipation.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[XD,YD] = meshgrid(xd,yd);
[V,gV] = Pot(XD,casV); 
pp=exp(-bet*(((YD.^2)/2)+V));
Dpp=(-bet)*(pp.*gV);
cn=dx*dy*sum(sum(pp(2:end,2:end)));                                        %% Integral of the theoretical invariant density on the cartesian product formed by xd and yd                                      
p=(1/cn)*pp;                                                               %% Normalization of the theoretical invariant density
Dp=(1/cn)*Dpp;                                                             %% derivative with respect to x of the normalized invariant density 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretization of the Drift of the SDE
[s] =S(XD,casS);                                                           %% compute of  the function s of the SDE
g=-((s.^2).*YD)-gV  ;                                                      %% preallocate for efficiency
DD=g.*p;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% discretization of the minus the gV  multiplied by the constant  -bet 
MbetgV=-bet*gV(1,:);                                                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute of the  Marginal invariant densities
px=dy*sum(p(2:end,:));                                                     %% theoretical marginal invariant density of X 
py=dx*sum(p(:,2:end),2);                                                   %% theoretical marginal invariant density of Y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

