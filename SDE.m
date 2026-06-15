function [UX,UY,hn]=SDE(n,gamm,rm,bet,x0,y0,casV,casS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates the approximation of a discretized trajectory of
% the process (X,Y)=(X_t,Y_t)_(t>=0)that resolves the damping hamiltonian 
% system defined by:  
%           dX_t=Y_tdt  
%           dY_t=sqrt(2/bet)*s(X_t)*dW_t -(s(X_t)^2*Y_t+ gV(X_t))*dt,
% where W is a Brownian motion and gV denotes the gradient of the potential
% V. The simulation will be denoted by (X,Y)_{k=1}^L. From the sample 
% (X,Y)_{k=1}^L is taken a subsample (UX,UY)_{i=1}^n and from this
% subsample will be computed the non parametric estimation of the invariant
% density of the process (X,Y).
%
% Problem parameters and input variables: n denotes the number of time 
% steps;rm determines the relation delta=hn/rm between the time steps delta 
% and hn; the constants bet and  gamm  are \beta and \gamma in the model; 
% the values x0 and y0 are the initial conditions of the sstem; casV and 
% casS determine the potential and the function s of the SDE that will be 
% used. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem parameters, input variables and output variables
L=rm*(n+1);                                                                %% L is the size of the  sample (X,Y)_{k=1}^L, rm should be entero.
X=zeros(1,L);                                                              %% preallocate for efficiency
Y=zeros(1,L);                                                              %% preallocate for efficiency
X(1)=x0 ;                                                                  %% x_0 is the initial condition for X
Y(1)=y0;                                                                   %% y_0 is the initial condition for Y
hn=n^(-gamm);                                                              %% hn is the time step corresponding to the sub-samples stored in (UX,UY)
delta=(1/rm)*hn;                                                           %% delta is the time step used to compute the sample (X,Y)_{k=1}^L. 
c=sqrt(2/bet);                                                             %% c is the constant that multiplies s in the diffusion coefficient of the system
%% Compute  X and Y. Lines of commands are constructed from the scheme in
%% page 96 of the book ''Free Energy Computations''                       
for j=1:L-1                                                            
    dW=sqrt(delta/2)*randn(1,2);                                           %% random numbers for the difussion term of the SDE
    % step 1                                                   
    [s] =S(X(j),casS);                                                     %% output of the function s
    cj14=1+(delta*(s^2)/4);                                                %% constant 
    Yj14=((1-(delta*(s^2)/4))*Y(j)+ (dW(1)*c*s))/cj14;                     %% approximate value of a trajectory of  Y in time t=j*h_n+1/4
    [ V,gV ] = Pot(X(j),casV);                                             %% output of the function Pot.m that compute the potential and its gradient
    H1=((1/2)*(Yj14^2))+ V;                                                %% value of the Hamiltonian in (X(i,j),Yj14)
    % step 2
    Yj12=Yj14-((delta/2)*gV);                                              %% approximate value of a trajectory of  Y in time t=j*h_n+1/2
    Xj1t=X(j)+(delta*Yj12);                                                %% temporal value of a trajectory of  X in time t=j*h_n+1
    [ V,gV ] = Pot(Xj1t,casV);                                             %% output of the function Pot.m that compute the potential and its gradient
    Yj34t=Yj12-((delta/2)*gV);                                             %% temporal value of a trajectory of  Y in time t=j*h_n+3/4
    % step 3
    H2=((1/2)*(Yj34t^2))+V;                                                %% value of the Hamiltonian in (Xj1t,Yj34t)
    r=exp(-bet*(H2-H1));                                                   %% r is the Metropolis-Hastings ratio
    mp=min(1,r);                                                           %% probability of acceptance or rejection
    u=rand;                                                                %% random value that determine the value of the trayectories in time t=(j+1)*h_n for X and in time t=h_n*(j+3/4)
    if (u<=mp)
       X(j+1)=Xj1t;                                                        %% possible  assignation of the value of X in the  time  t=(j+1)*h_n
       Yj34=Yj34t;                                                         %% possible  assignation of the value of Y in the fractional  time t=jh_n+(3/4)
    else
       X(j+1)=X(j);                                                        %% possible  assignation of the value of X in the  time  t=(j+1)*h_n
       Yj34=-Yj14;                                                         %% possible  assignation of the value of Y in the fractional  time t=j*h_n+(3/4)
    end
    %step 4                                                       
    [s] =S(X(j+1),casS);
    cj1=1+(delta*(s^2)/4); 
    Y(j+1)=((1-(delta*(s^2)/4))*Yj34+(c*s*dW(2)))/cj1;   %%%%              %% assignation of the value of Y in the  time  t=(j+1)h_n
end

%% Compute  (UX,UY)
UX=X(rm:rm:L);                                                             %% UX is a sub-samplig of X
UY=Y(rm:rm:L);                                                             %% UY is a sub-sampling of Y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end