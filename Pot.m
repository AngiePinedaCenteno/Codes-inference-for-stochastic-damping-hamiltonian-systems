function [ V,gV ] =Pot(position,casV)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Potential_Hypoelyptic_Diffusion=casV;
switch Potential_Hypoelyptic_Diffusion
    case 1 
        %% Potential 1 (Trigonometric) 
        V=-cos(position);                                                  %% values of the potential on position
        gV=sin(position);                                                  %% values of the gradient of the potential
    case 2 
        %% Potential 2 (Trigonometric) 
        V=2-cos(position);                                                 %% values of the potential on position
        gV=sin(position);                                                  %% values of the gradient of the potential
    case 3 
        %% Potential 3 (Polynomial )
         V=position.^2;
         gV=2*position;
    case 4
        %%  Potential 4 (Polynomial)
        V=(position.^4)-(2*(position.^2))+4;
        gV=(4*(position.^3))-(4*position);
    otherwise
        %%  potential  5 (Polynomial)
        V=(6*(position.^6))-(12*(position.^4))+(6*(position.^2))+(1/7);
        gV=(36*(position.^5))-(48*(position.^3))+(12*position);
end
end