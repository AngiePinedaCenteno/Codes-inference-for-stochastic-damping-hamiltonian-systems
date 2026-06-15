function [driftE0And1]=DensityAndDrift_0And1(UX,UY,xd,sxd,n,b1n,b2n,hn,varep) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the estimate of the invariant density and the drift 
% on the ordered pairs (x,1) and (x,0), for x in a discretized domain deter-
% mined by the sample which approximate a trajectory of the position.
% The input variables of this function are described in the function
% Inputs_FluctuationDissipation.m 
% The outputs of this function are decribes in the function
% Outputs_FluctuationDissipation.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute  of the estimation of the invariant density
yd=[0,1];                                                                  %% discretization of the domain of the dendity in y -direction
syd=length(yd);                                                            %% size of discretization of the domain of the density in the  y -direction
driftE0And1=zeros(syd,sxd);                                                %% preallocate for efficiency
cx=1/b1n; cy=1/b2n;
Kpx=zeros(1,n);                                                            %% preallocate for efficiency
Kpy=zeros(syd,n);                                                          %% preallocate for efficiency
DUY=(1/hn)*(UY(2:end)-UY(1:n));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
for i=1:syd               
    YI=cy*(yd(i)-UY(1:n));                                                 %% argument of the kernel in y- direction
    for l=1:n
        [Key] = EK(YI(l));
        Kpy(i,l)=Key;                                                      %% value of the Epanechnikov kernel associated to the point yd(i)
    end
end
%%%%% Estimation of the invariant density,derivative anddrift %%%%%%%%%%%%%
for j=1:sxd
    XI=cx*(xd(j)-UX(1:n));                                                 %% argument of the kernel in X
    for l=1:n
        [Kex] = EK(XI(l));
        Kpx(l)=Kex;                                                        %% value of the Epanechnikov kernel associated to the point xd(j)
    end    
    for i=1:syd
    K=Kpx.*Kpy(i,:);
    Kh=K.*DUY;
    dKp=(1/(n*b1n*b2n))*sum(K);                                            %% value of the estimate of the invariant density in the point (xd(j),yd(i))
    Hn=(1/((n-1)*b1n*b2n))*sum(Kh(1:n-1));
    driftE0And1(i,j)=Hn/(dKp+varep);                                       %% value of the estimate of the drift in the point (xd(j),yd(i))
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end