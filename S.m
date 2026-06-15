function [s] =S(position,casS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the values of the function s of the SDEdefined in
% SDE.m.
%
% Input variables: position is the element where will be evaluated the
% s function, and casS is the variable which determine the s function that
% will be used to compute the the approximate solution of the SDE defined
% in SDE.m.
% Ouput variables: s is the value of the function S evaluated in the
% value,vector or matrix position.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_Hypoelyptic_Diffusion=casS;
switch S_Hypoelyptic_Diffusion
    case 1 
        [ei,ed]=size(position);
        s=2*ones(ei,ed);                                                   %% sigma constant
    case 2
        s=exp(-1./((position.^2)+1));                                      %% sigma Not constant
    case 3   
        s=10*exp(-1./((position.^2)+1));                                      
    case 4
        s=sqrt(2+cos(5*position));                                         
    case 5
        s=sqrt(2+cos(8*position));                                          
    case 6
        s=sqrt(2+cos(5*position)); 
    case 7   
        s=exp(2)*exp(-1./((position.^2)+1));
    case 8   
        s=sqrt(2+cos(3*position)); 
    case 9   
        s=sqrt(2+cos(10*position));
    case 10                                                                 
        s=5*sqrt(2+cos(10*position));  
    case 11                                                                 
        s=sqrt((1+0.0001)+cos(10*position)); 
    case 12                                                                 
        s=sqrt((1+0.0001)+cos(position)); 
    case 13                                                                 
        s=exp(-1./((position.^2)+3));  
     case 14                                                                 
        s=15*exp(-1./((position.^2)+1));  
     case 15                                                                 
        s=5*exp(-1./((position.^2)+1));  
     case 16                                                                 
        s=3*exp(-1./((position.^2)+1));  
     case 17                                                                 
        s=(1.5)*exp(-1./((position.^2)+1));  
    case 18                                                                 
        s=(0.5)*exp(-1./((position.^2)+1));  
    case 19
        s=sqrt(2+cos(position)); 
    otherwise
        s=sqrt(1+cos(4*position));                                         
end
end