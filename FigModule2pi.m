function FigModule2pi(n,gamm,varep,xd,yd,alpha1,alpha2,p,px,py,dKp,dKpx,dKpy,Dp,DdKp,g,DD,driftE,driftE0And1,Hn,EbetgVM,MbetgV,casV,casS,casGS) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function FiguresHypoelypticDiffusionFull.m generates the graphs of:
% Case 1: a discretized invariant density and its estimated version, the
%         discretizations of the marginal invariant densities and  their 
%         estimations, a discretization of the Derivative with respect to x
%         of the Invariant Density  and its estimate, a discretization of
%         the theoretical drift and its estimate and a discretization of 
%         the function $H_n$ that approximate the drift multiplied by the
%         invariant density.These graphics are obtained when casGS, in the 
%         input parameters, is equal to 1.
% Case 2: the estimate of $-\beta\nablaV(x)$ for fixed values of y. The 
%         values of y used have the condition |y|<=1.7 when the potential 
%         is polynomial and |y|<=1.5 when the potential is periodic .For 
%         these values of y the estimate of the marginal invariant density 
%         for the velocity is not zero.
%         Also,for the graphics there  a restriction for the values of x. 
%         The variable x is in a rectangular discretization of [-0.8,0.8] 
%         when the potential is polynomial, and the variable x is in a rec-
%         tangular discretization of [-2.5,2.5] when the potential is pe-
%         riodic. These restriction for the values of x and y are used for
%         the representations in the rest of the cases
%         These graphics are obtained when casGS,in the input parameters,is
%         equal to 2.
% Case 3: the projections of the drift for fixed values of y and its esti-
%         mate. These graphics are obtained when casGS in the input parame-
%         ters is equal to 3.
% Case 4: the projections of the drift for fixed values of x and its 
%         estimate.
%         These graphics are obtained when casGS in the input parameters is
%         equal to 4.
% Case 5: the comparative graphs of the theoretical variance versus their 
%         estimate. The first figure that you can see correspond to the 
%         comparison between the variance and its estimation for y=1.After,
%         you must press the enter key to obtain the comparison between the
%         variance and its estimation from several  fixed values of y.These
%         graphics are obtained when casGS in the input paameters is equal 
%         to 5.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ((casV==1) || (casV==2) )
    AD=find(abs(xd)<=2.5);
    ADy=find(abs(yd)<=1.5);
else
    AD=find(abs(xd)<=0.8);
    ADy=find(abs(yd)<=1.7);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figures
%%%%%%%%%%%%%%%%%%%% Estimation  Invariant Density %%%%%%%%%%%%%%%%%%%%%%%%
GS_Hypoelyptic_Diffusion=casGS;
switch GS_Hypoelyptic_Diffusion
    case 1
        close all
        %ADy=find(abs(yd)<=1.5);
        %%%%%%%%%%%%%% Drift and its estimate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%  Drift       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        li1=min(min(g(ADy,AD)));li2=min(min(driftE(ADy,AD)));li=min(li1,li2);
        ls1=max(max(g(ADy,AD)));ls2=max(max(driftE(ADy,AD)));ls=max(ls1,ls2);        
        figure(5)
        subplot(1,2,1);
        zlim([li ls]);
        set(gca,'NextPlot','replaceChildren');
        surfl(xd(AD),yd(ADy),g(ADy,AD));                                                            
        shading interp
        colormap gray
        title(['Drift Adjusted Domain :' '   gamm= ',num2str(gamm) , ', alpha1= ', num2str(alpha1), ', alpha2= ', num2str(alpha2) ]);
        xlabel('position')
        ylabel('velocity')
        view([60, 70, 20]);
        %%%%%%%%%%%%%%%%%% Estimation of the Drift Adjusted Domain %%%%%%%%
        subplot(1,2,2);
        zlim([li ls]);
        set(gca,'NextPlot','replaceChildren');
        surfl(xd(AD),yd(ADy),driftE(ADy,AD));                                                                                                                   
        shading interp
        colormap gray
        title(['Estimated Drift  Adjusted Domain :' '   gamm= ',num2str(gamm) , ', alpha1= ', num2str(alpha1), ', alpha2= ', num2str(alpha2) ]);
        xlabel('position')
        ylabel('velocity')
        view([60, 70, 20]);
        %%%%%%%%%%%%% Drift multiply by density  and Hn %%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%% Drift multiply by density          %%%%%%%%%%%%%%%%%%
        li1=min(min(DD));li2=min(min(Hn));li=min(li1,li2);
        ls1=max(max(DD));ls2=max(max(Hn));ls=max(ls1,ls2);
        figure(4)
        subplot(1,2,1);
        zlim([li ls]);
        set(gca,'NextPlot','replaceChildren');
        surfl(xd,yd,DD);                                                            
        shading interp
        colormap gray
        title(['Drift multiply by density :' '   gamm= ',num2str(gamm) , ', alpha1= ', num2str(alpha1), ', alpha2= ', num2str(alpha2) ]);
        xlabel('position')
        ylabel('velocity')
        view([60, 70, 20]);
        %%%%%%%%%%%%%%%%%% Hn %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(1,2,2);
        zlim([li ls]);
        set(gca,'NextPlot','replaceChildren');
        surfl(xd,yd,Hn);                                                                                                                   
        shading interp
        colormap gray
        title(['Hn :' '   gamm= ',num2str(gamm) , ', alpha1= ', num2str(alpha1), ', alpha2= ', num2str(alpha2) ]);
        xlabel('position')
        ylabel('velocity')
        view([60, 70, 20]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Discretized derivative with respect to x of the invariant density  
        li1=min(min(Dp));li2=min(min(DdKp));li=min(li1,li2);
        ls1=max(max(Dp));ls2=max(max(DdKp));ls=max(ls1,ls2);
        figure(3)
        subplot(1,2,1);
        zlim([li ls]);
        set(gca,'NextPlot','replaceChildren');
        surfl(xd,yd,Dp);                                                                                                                   
        shading interp
        colormap gray
        title(['Derivative of the Invariant Density with respect to x :' '   gamm= ',num2str(gamm) , ', alpha1= ', num2str(alpha1), ', alpha2= ', num2str(alpha2) ]);
        xlabel('position')
        ylabel('velocity')
        %%%Estimation with respect to x of the derivative of the invariant 
        %%% density 
        subplot(1,2,2);
        zlim([li ls]);
        set(gca,'NextPlot','replaceChildren');
        surfl(xd,yd,DdKp);                                                                                                                   
        shading interp
        colormap gray
        title(['Estimate of the Derivative of the Invariant Density with respect to x :' '   gamm= ',num2str(gamm) , ', alpha1= ', num2str(alpha1), ', alpha2= ', num2str(alpha2) ]); 
        xlabel('position')
        ylabel('velocity')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%% Marginal Density with respect to X  %%%%%%%%%%%%%
        figure(2)
        subplot(1,2,1);
        plot(xd,px,'r'),hold on
        plot(xd,dKpx,'b--'), hold off
        title(['Invariant Density for X:' '   gamm= ',num2str(gamm) , ', alpha1= ', num2str(alpha1), ', alpha2= ', num2str(alpha2) ]); 
        %%%%%%%%%%%%%%%%% Marginal Density with respect to Y  %%%%%%%%%%%%%
        subplot(1,2,2);
        plot(yd,py,'r'),hold on
        plot(yd,dKpy,'b--'), hold off
        title(['Invariant Density for Y:' '  gamm= ',num2str(gamm) , ', alpha1= ', num2str(alpha1), ', alpha2= ', num2str(alpha2) ]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% Discretized Invariant Density %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        li1=min(min(p));li2=min(min(dKp));li=min(li1,li2);
        ls1=max(max(p));ls2=max(max(dKp));ls=max(ls1,ls2);
        figure(1)
        subplot(1,2,1);
        zlim([li ls]);
        set(gca,'NextPlot','replaceChildren');
        surfl(xd,yd,p);                                                            
        shading interp
        colormap gray
        title(['Invariant Density:' '  n= ',num2str(n), ', gamm= ',num2str(gamm) , ', alpha1= ', num2str(alpha1), ', alpha2= ', num2str(alpha2) ]); 
        xlabel('position')
        ylabel('velocity')
        subplot(1,2,2);
        zlim([li ls]);
        set(gca,'NextPlot','replaceChildren');
        surfl(xd,yd,dKp);                                                                                                                   
        shading interp
        colormap gray
        title(['Estimation Invariant Density:' '  gamm= ',num2str(gamm) , ', alpha1= ', num2str(alpha1), ', alpha2= ', num2str(alpha2) ]);
        xlabel('position')
        ylabel('velocity')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% The estimates of $-\beta\nabla V(x)$ for fixed values of y %%%%%%%
    case 2
        %%%%%%%%%%%  Adjusted domain  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        close all
        %ADy=find(abs(yd)<=1.1);
        li=min(min(EbetgVM(:,AD)));
        ls=max(max(EbetgVM(:,AD)));
        for aa=ADy %%%for aa=1:lADy
            y=yd(aa); %%%y=yd(ADyM(aa));
            figure(aa)
            ylim([li ls]);
            set(gca,'NextPlot','replaceChildren');
            plot(xd(AD),EbetgVM(aa,AD),'b--')                                                  
            hold on
            plot(xd(AD),MbetgV(AD),'r')
            hold off
            title(['Minus gradient of the Potential and its estimate on an Adjusted Domain for' '  y= ',num2str(y) ]);
            xlabel(['Parameters:' ' gamm= ',num2str(gamm), ', alpha1= ',num2str(alpha1) , ', alpha2= ', num2str(alpha2), ', n= ', num2str(n) ])
            pause(1)
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%  Drift for fixed values of y   %%%%%%%%%%%%%%%%%%%%%%
    case 3
        close all
        ADy=find(abs(yd)<=1.1);
        li=min(min(g(ADy,AD)));                                            
        ls=max(max(g(ADy,AD)));                                            
        for aa= ADy
            y=yd(aa);
            figure(aa)
            ylim([li ls]);
            set(gca,'NextPlot','replaceChildren');
            plot(xd(AD),driftE(aa,AD),'b--')                                                  
            hold on
            plot(xd(AD),g(aa,AD),'r')
            hold off
            %title(['Drift and its estimate on an adjusted domain for' '  y= ',num2str(y)  ]);
            title(['Drift and its estimate for the fixed value' '  y= ',num2str(y) 'on an adjusted domain' ]);
            xlabel(['Parameters:' ' gamm= ',num2str(gamm), ', alpha1= ',num2str(alpha1) , ', alpha2= ', num2str(alpha2), ', n= ', num2str(n) ])
            pause(1)
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% Drift for fixed values of x %%%%%%%%%%%%%%%%%%%%%%%%
      case 4
        close all
        ADy=find(abs(yd)<=1.5);
    %%%%%%%%%%%  Adjusted domain  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        li=min(min(g(ADy,AD)));
        ls=max(max(g(ADy,AD)));
        for aa=AD
            x=xd(aa);
            figure(aa)
            ylim([li ls]);
            set(gca,'NextPlot','replaceChildren');
            plot(yd(ADy),driftE(ADy,aa)','b--')                                                  
            hold on
            plot(yd(ADy),g(ADy,aa),'r')
            hold off
            title(['Drift and  its estimate for a fixed value of ' '  x= ',num2str(x) 'on an Adjusted domain']);
            xlabel(['Parameters:' ' gamm= ',num2str(gamm), ', alpha1= ',num2str(alpha1) , ', alpha2= ', num2str(alpha2), ', n= ', num2str(n) ])
            pause(1)
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Estimation of the variance from Two Independent samples %%%%%%
    case 5
        close all
        se=-(driftE0And1(2,:)- driftE0And1(1,:));
        [s] =S(xd,casS);
        sets=s.^2;
        figure(1)
        ylim([-15 15]);
        set(gca,'NextPlot','replaceChildren');        
        plot(xd,se,'b--')
        hold on
        plot(xd,sets,'r')
        hold off                                                           
        title('Estimation of s(x)*s(x)^T ' );
        xlabel(['Parameters:' '  varep= ', num2str(varep) , ', n= ', num2str(n)])
        pause
    %%%%%%%%%%%%%% Estimation of -s(x)*s(x)^T*y   %%%%%%%%%%%%%%%%%%%%%%%%%
        close all
        [s] =S(xd,casS);
        sets=s.^2;
        %ADy=find(abs(yd)<=1.1);
        for aa=ADy
            y=yd(aa);
            se=-(driftE(aa,:)- driftE0And1(1,:))/y;
            figure(1+aa)
            ylim([-15 15]);
            set(gca,'NextPlot','replaceChildren');
            plot(xd(AD),se(AD),'b--')
            hold on
            plot(xd(AD),sets(AD),'r')                                                          
            title([ '  y= ', num2str(y) ]);
            xlabel(['Parameters:' '  varep= ', num2str(varep) , ', n= ', num2str(n)])
            hold off 
            pause(0.5)
        end
end
end

