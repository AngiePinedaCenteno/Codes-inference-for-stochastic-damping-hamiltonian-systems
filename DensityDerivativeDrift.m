function [dKp,dKpx,dKpy,driftE,DdKp,Hn,xd,yd,sxd,syd]=DensityDerivativeDrift(UX,UY,b1n,b2n,hn,dx,dy,n,varep)   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes non parametric estimates of the invariant density 
% of the process (X,Y)=(X_t,Y_t)_(t>=0) which resolves the system:   
%                  dX_t=Y_tdt  
%                  dY_t=sqrt(2/bet)*s*dW_t -(s^2*Y_t+gV(X_t))*dt,
% where: bet is a constant, W is a Brownian motion and gV denotes the 
% gradient of the potencial V.
% Also,this function computes non parametric estimations of the drift term, 
% and the derivative with respect to x of the invariant density.
% The input variables of this function are described in the function
% Inputs_FluctuationDissipation.m 
% The outputs of this function are decribes in the function
% Outputs_FluctuationDissipation.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretized domain of the invariant density, derivative and drift term
xmind=min(UX(1:n)-b1n);                                                    %% minimum value of x for the domain                                        
xmaxd=max(UX(1:n)+b1n);                                                    %% maximum values of x for the domain
ymind=min(UY(1:n)-b2n);                                                    %% minimum value of y for the domain
ymaxd=max(UY(1:n)+b2n);                                                    %% maximum value of y for the domain
xd=xmind:dx:xmaxd;                                                         %% discretization of the domain of the dendity in the x -direction
yd=ymind:dy:ymaxd;                                                         %% discretization of the domain of the dendity in y -direction
sxd=length(xd);                                                            %% size of discretization of the domain of the dendity in the x -direction
syd=length(yd);                                                            %% size of discretization of the domain of the dendity in the  y -direction
dKp=zeros(syd,sxd);                                                        %% preallocate for efficiency                                                   
dKpx=zeros(1,sxd);                                                         %% preallocate for efficiency
Hn=zeros(syd,sxd);                                                         %% preallocate for efficiency
driftE=zeros(syd,sxd);                                                     %% preallocate for efficiency
DdKp=zeros(syd,sxd);                                                       %% preallocate for efficiency
cx=1/b1n; cy=1/b2n;                                                        %% coeficients of the terms of the estimator of the density
Kpx=zeros(1,n);                                                            %% preallocate for efficiency
Kpy=zeros(syd,n);                                                          %% preallocate for efficiency
KGpy=zeros(syd,n);                                                         %% preallocate for efficiency
DUY=(1/hn)*(UY(2:end)-UY(1:n));                                            %% Factors of terms in the estimator H_n
%%%%%%% Estimation of the marginal invariant density for the velocity %%%%%    
for i=1:syd               
    YI=cy*(yd(i)-UY(1:n));                                                 %% arguments of the kernel in the density estimator
    for l=1:n
        [Key] = EK(YI(l));
        Kpy(i,l)=Key;                                                      %% value of the Epanechnikov kernel associated to the point yd(i)
    end
    KGpy(i,:)=normpdf(YI,0,1);                                             %% values of the Gaussian kernel associated to the point yd(i)                    
end
dKpy=(cy/n)*sum(Kpy,2);                                                    %% estimate of the marginal invariant density for  the velocity 
% Estimation of the invariant density, marginal invariant density for the 
% position, derivative and drift 
for j=1:sxd
     XI=cx*(xd(j)-UX(1:n));                                                %% argument of the kernel in X
    for l=1:n
        [Kex] = EK(XI(l));
        Kpx(l)=Kex;                                                        %% value of the Epanechnikov kernel associated to the point xd(j)
    end
    dKpx(j)=(cx/n)*sum(Kpx);                                               %% value of the  estimate of the marginal invariant density for velocity  in the point xd(j)
    for i=1:syd
    K=Kpx.*Kpy(i,:);                                                       %% terms of the estimate of the invariant density  for the point (xd(j),yd(i))
    Kh=K.*DUY;                                                             %% terms of the drift estimate for the point (xd(j),yd(i))
    Kd=(-1)*((KGpy(i,:).*normpdf(XI,0,1)).*XI);                            %% terms of the estimate of the derivative with respect to x of the invariant density in the point (xd(j),yd(i))
    dKp(i,j)=(1/(n*b1n*b2n))*sum(K);                                       %% value of the estimate of the invariant density in the point (xd(j),yd(i))
    Hn(i,j)=(1/((n-1)*b1n*b2n))*sum(Kh(1:n-1));
    driftE(i,j)=Hn(i,j)/(dKp(i,j)+varep);
    DdKp(i,j)=(1/(n*(b1n^2)*b2n))*sum(Kd);                                 %% values of the estimator of the derivative with respect to x of the invariant density in the point (xd(j),yd(i))
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end