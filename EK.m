function [Ke] = EK(w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function compute the value 'Ke' of the Epanechnikov kernel in the 
% point w
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (abs(w)<= 1) 
        Ke=(3/4)*(1-(w^2));                                                %% possible value of the Epanechnikov kernel in the point w 
    else 
        Ke=0;                                                              %% possible value of the Epanechnikov kernel in the point w
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

