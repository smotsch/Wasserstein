% This code computes the distance between the particle-level equation
% and the macro-level equation. Author: Shane Lubold (shane.lubold@asu.edu)
% Vel1 corresponds to rhoMacro and is a deterministic quantity. Vel2
% corresponds to rhoMicro and is a random quantity.

clc; clear; close all;

% DeltaX = [2 1 .5 .25 .0125];
% BigT = [10 10^2 10^3 10^4 10^5]
% for m = 1:length(DeltaX)


num = [5 10 10^2 10^3 10^4]; 
dx = [2 1 .5 .25 .125]; % number particles


for i = 1:length(num)
    a   = -5;                              % lower bound on grid
    b   = 5; 
    dx1 = dx(i);
    intX = a:dx1:b; 
    IC = normpdf(intX,0,1) ;
    D  = 1;  
    dt = .01;
    T =  0:dt:1;  
    rhoMacro  = KPP_Macro(IC,intX,D,dt,T,2);
    for k = 1:1000
        %% Parameter model
        D  = 1;                                % noive level
        sigma = sqrt(2*D);
        l0  = .2;                                % radius Initial condition
        a   = -5;                              % lower bound on grid
        b   = 5;                               % upper bound on grid
        dx1 = dx(i);                                % meshsize for micro
        dx3 = .1;                             % meshsize for Hist for rhoMicro
        dx2 =dx1;                                % meshsize for macro
        intX = a:dx1:b;                        % physical discretization
        dt = .01;                              % delta t
        T =  0:dt:1;                           % time vector
        IC = normpdf(intX,0,1) ;                               %  Initial condition 

        %% Solve  Macro
        % rhoMacro1 = MacroFinal(intX,dx2,T,dt,sigma,a,b,l0);


        %% Solve  Micro

        xLoc1(1) = 0; xLoc2(1) = 0;             % Initialize position vectors


        %[rhoMicro,FINALX] = MicroFinalFast(IC,intX,D,dt,T,2);
        [rhoMicro,X] = KPP_Micro(IC,intX,D,dt,T,num(i));
        tic
        %       for j = 2:length(T)
        %           plot(intX,rhoMacro(:,j),'ko',intX,rhoMicro(:,j),'g'); 
        %           legend('Macro','Micro'); grid on;   axis([a b 0 1]);
        %           pause(.01);
        
        WD_cont(i,k) = WD_Cont(rhoMicro(:,end),rhoMacro(:,end),intX);
        %               WD_discreet(j)= WD_discreet_cont(X,intX,rhoMacro(:,end));
        %              xLoc1(j) = CalculateVelocity(rhoMicro(:,j),intX,xLoc1(j-1));
        %              plot(abs(FINV(:,j)-GINV(:,j))); grid on; hold on;
        %       end


        

        close all;
    end
    % toc
end


% figure; plot(Vel1(2:end),'k--'); grid on; 
% figure; plot(Vel2(2:end),'r--'); grid on; 

%How does WD({x_i},hist) change with \Delta x.
% 
% [WDStepCompare(k), WDHistCompare(k)] = CombinedWDCalculator(FINALX{j},intX,rhoMicro(:,1),rhoMicro(:,1),dx1);
% end
% WD_contFinal(m) = (WD_cont(end));
% LogWD_contFinal(m) = log(WD_cont(end));
% end
% 
% % fun = @(x,xdata) x(1).*(1./xdata).^(x(2));
% % x0 = [.3,1];
% % times = linspace(BigT(1),BigT(end));
% % fit = lsqcurvefit(fun,x0,BigT,WD_contFinal);
% 
% fit = polyfit(DeltaX,WD_contFinal,1);
% yfit = fit(2) + fit(1)*DeltaX;
% 
% 
% plot(DeltaX,WD_contFinal,'ko'); hold on; plot(DeltaX,yfit,'k-')
% leg = legend('Data','Fitted Model'); grid on;
% set(leg,'position',[.2 .57 .3 .3])
% xlabel('\Delta x'); ylabel('WD'); title('WD vs. \Delta x');
% 
% 
% 
% Fit = polyfit(BigT,WD_contFinal,2);
% Fit = polyval(Fit,BigT);
% FitLog = polyfit(BigT,LogWD_contFinal,1);
% figure; plot(BigT,WD_contFinal); grid on;
% hold on; plot(BigT,Fit); hold off;
% xlabel('N'); ylabel('WD');


