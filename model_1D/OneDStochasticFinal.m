% This code computes the distance between the particle-level equation
% and the macro-level equation. Author: Shane Lubold (shane.lubold@asu.edu)
% Vel1 corresponds to rhoMacro and is a deterministic quantity. Vel2
% corresponds to rhoMicro and is a random quantity.

clc; clear; close all;

%% Parameter model
D  = 1;                                % noive level
sigma = sqrt(2*D);
l0  = .2;                                % radius Initial condition
num = [10^4];                          % number particles
a   = -5;                              % lower bound on grid
b   = 5;                               % upper bound on grid
dx1 =.05;                                % meshsize for micro
dx3 = .1;                             % meshsize for Hist for rhoMicro
dx2 =.05;                                % meshsize for macro
intX = a:dx1:b;                        % physical discretization
dt = .001;                              % delta t
T =  0:dt:1;                           % time vector
IC = normpdf(intX,0,1) ;                               %  Initial condition 

%% Solve  Macro
% rhoMacro1 = MacroFinal(intX,dx2,T,dt,sigma,a,b,l0);
tic
rhoMacro  = KPP_Macro(IC,intX,D,dt,T,2);
toc
%% Solve  Micro
tic
xLoc1(1) = 0; xLoc2(1) = 0;             % Initialize position vectors

for i = 1:length(num)
    %[rhoMicro,FINALX] = MicroFinalFast(IC,intX,D,dt,T,2);
    [rhoMicro,X] = KPP_Micro(IC,intX,D,dt,T,num(i));
    tic
      for j = 2:length(T)
%           plot(intX,rhoMacro(:,j),'ko',intX,rhoMicro(:,j),'g'); 
%           legend('Macro','Micro'); grid on;   axis([a b 0 1]);
%           pause(.01);
        

%             [HistCDF,intXHist] = cdf_clean(intX,rhoMicro(:,j)');
%             [MacroCDF,intXMacro] = cdf_clean(intX,rhoMacro(:,j)');
             [WD_cont(j),MacroCDF(:,j),xMac] = WD_Cont(rhoMicro(:,j),rhoMacro(:,j),intX);
             [WD_discreet(j),F(:,j),FINV(:,j),GINV(:,j)] = WD_discreet_cont(X,MacroCDF(:,j),xMac);
%             [WDHist(j)] = WDFinal(HistCDF,MacroCDF,intXHist,intXMacro);
%             [WDParticle(j),pCombined,FINV,GINV] = WDFinal(StepCDF,MacroCDF,FINALX{j},intXMacro);
%            [Vel1(i,j), Vel2(i,j),xLoc1(j),xLoc2(j)] = CalculateVelocity(T,rhoMacro(:,j),rhoMicro(:,j),intX,xLoc1(j-1),xLoc2(j-1));
%             plot(abs(FINV(:,j)-GINV(:,j))); grid on; hold on;
      end
     toc
     hold off;
end
toc

% figure; plot(Vel1(2:end),'k--'); grid on; 
% figure; plot(Vel2(2:end),'r--'); grid on; 

%How does WD({x_i},hist) change with \Delta x.
% 
% [WDStepCompare(k), WDHistCompare(k)] = CombinedWDCalculator(FINALX{j},intX,rhoMicro(:,1),rhoMicro(:,1),dx1);
% end



