%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%              - Exercício Computacional de MP208 -             %
%    --- Optimal Filtering with Aerospace Applications ---      %
%                                                               %
%              Autor: João Filipe R. P de A. Silva              %
%                                                               %
%              Main Script: Análise de Estimadores              %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cache Clean-up
clear all
close all
clc

%% General Variables declaration

est.R = 0.01;               %Measurement Noise variance
est.Q = diag([0.01 0.04]);  %States Disturbance variance
est.P = diag([1E-4 1E-8]);  %X Covariance Matrix
est.m = [1 0]';             %X Expected Value
est.e1 = [1 0]';            %Unitary Vector
est.e2 = [0 1]';            %Unitary Vector
est.T = 0.1;                %Sampling Time
est.x = [est.m];            %Estimatior initial Value
t = 0:est.T:20;             %Simulation Period
N = 10;                     %Number of Realizations in the simulation
est.mean = zeros(2,length(t));          %Error Sample Mean (Kalman Filter)
est.RMS = zeros(2,length(t));           %RMS Error (Kalman Filter)
est.Infmean = zeros(2,length(t));       %Error Sample Mean (Information Filter)
est.InfRMS = zeros(2,length(t));        %RMS Error (Information Filter)
Pkk = zeros(2,2,N);             %Estimation Error Covariance Matrix for each realization

%% Plant Modeling

sys.A = [1 0.1; 0 1];                      %System States Matrix
sys.B = [0.005; 0.1];                      %System Input Matrix
sys.C = [1 0];                             %System Output Matrix
sys.y = 0;                                 %System Output Initialization;
sys.y_bar = 5;                             %Command Input Initialization;

sys.P = est.P;                             
sys.m = est.m; 

%% Simulation

for k = 1:N
    
    sys.x = sqrtm(sys.P)*randn(2,1) + sys.m;             %Initial value of X
    est.P = diag([1E-4 1E-8]);                           %Initial Estimate Covariance Matrix
    est.x = [est.m];                                     %Initial Estimate of X
    
    for cont = 1:length(t)

        %Model Behaviour

        sys.u = 10*(sys.y_bar - est.e1'*sys.x) - 2*est.e2'*sys.x;   %Control Input
        sys.w = sqrtm(est.Q)*randn(2,1);                %State Disturbance Realization
        sys.v = sqrt(est.R)*randn(1,1);                 %Measurement Noise Realization

        sys.x = sys.A*sys.x + sys.B*sys.u + sys.w;      %System states dynamic
        sys.y = sys.C*sys.x + sys.v;                    %System measurement dynamic

        y(k,cont) = sys.y;                              %Output plotting variable;
        y_bar(k,cont) = sys.y_bar;                      %Command Input plotting variable;
        
        % Information Filter

        est = InfFil(sys,est);                          %Information Filter function
        exInf = inv(est.L)*est.z - sys.x;               %Information Filter estimation Error

        exInf1(k,cont) = est.e1'*exInf;                 %Isolating exInf component 1
        exInf2(k,cont) = est.e2'*exInf;                 %Isolating exInf component 2
        
        est.Infmean(1,cont) = est.Infmean(1,cont) + exInf(1);          %Sample mean of IF Estimation Error of x1
        est.Infmean(2,cont) = est.Infmean(1,cont) + exInf(2);          %Sample mean of IF Estimation Error of x2
         
        est.InfRMS(1,cont) = est.InfRMS(1,cont) + exInf(1)*exInf(1)';     %Sum of Squared Errors in exInf(1)
        est.InfRMS(2,cont) = est.InfRMS(2,cont) + exInf(2)*exInf(2)';     %Sum of Squared Errors in exInf(2)
        
        est.InfPplot(1,cont) = sqrtm(inv(est.L(1,1)));      %Theoretical IF Standard Deviation of exInf(1)
        est.InfPplot(2,cont) = sqrtm(inv(est.L(2,2)));      %Theoretical IF Standard Deviation of exInf(2)
        
        % Kalman Filter

        est = KFil(sys,est);                            %Kalman Filter function
        ex = est.x - sys.x;                             %Kalman Filter estimation Error
        
        ex1(k,cont) = est.e1'*ex;                       %Isolating ex component 1
        ex2(k,cont) = est.e2'*ex;                       %Isolating ex component 1
        
        est.mean(1,cont) = est.mean(1,cont) + ex(1);          %Sample mean of KF Estimation Error of x1
        est.mean(2,cont) = est.mean(1,cont) + ex(2);          %Sample mean of KF Estimation Error of x2

        est.RMS(1,cont) = est.RMS(1,cont) + ex(1)*ex(1)';     %Sum of Squared Errors in ex(1)
        est.RMS(2,cont) = est.RMS(2,cont) + ex(2)*ex(2)';     %Sum of Squared Errors in ex(2)
        
        est.Pplot(1,cont) = sqrtm(est.P(1,1));          %Theoretical KF Standard Deviation of ex(1)
        est.Pplot(2,cont) = sqrtm(est.P(2,2));          %Theoretical KF Standard Deviation of ex(2)
         
%         

    end
    figure(1);                            
    plot(t,y(k,:));                     
    hold on;
    title('\textbf{Measured Outputs $y_k$ vs Command Input $\bar{y}_k$}','interpreter','LaTeX');
    xlabel('\textbf{Time [s]}','interpreter','LaTeX');
    ylabel('\textbf{Amplitude [\,-\,]}','interpreter','LaTeX');
    
    figure(2)
    plot(t,ex1(k,:));
    hold on;
    title('\textbf{Kalman Filter Estimation Error of state $x_1$}','interpreter','LaTeX');
    xlabel('\textbf{Time [s]}','interpreter','LaTeX');
    ylabel('\textbf{Amplitude [\,-\,]}','interpreter','LaTeX');
    
    figure(3)
    plot(t,ex2(k,:));
    hold on;
    title('\textbf{Kalman Filter Estimation Error of state $x_2$}','interpreter','LaTeX');
    xlabel('\textbf{Time [s]}','interpreter','LaTeX');
    ylabel('\textbf{Amplitude [\,-\,]}','interpreter','LaTeX');
    
    figure(4)
    plot(t,exInf1(k,:));
    hold on;
    title('\textbf{Information Filter Estimation Error of state $x_1$}','interpreter','LaTeX');
    xlabel('\textbf{Time [s]}','interpreter','LaTeX');
    ylabel('\textbf{Amplitude [\,-\,]}','interpreter','LaTeX');
    
    figure(5)
    plot(t,exInf2(k,:));
    hold on;
    title('\textbf{Information Filter Estimation Error of state $x_2$}','interpreter','LaTeX');
    xlabel('\textbf{Time [s]}','interpreter','LaTeX');
    ylabel('\textbf{Amplitude [\,-\,]}','interpreter','LaTeX');
    

end

est.mean(1,:) = est.mean(1,:)/N;                        %KF Estimation Error Sample Mean of x1
est.mean(2,:) = est.mean(2,:)/N;                        %KF Estimation Error Sample Mean of x2

est.pRMS(1,:) = sqrt(est.RMS(1,:)/(N-1));                 %KF Root Mean Square Estimation Error of x1
est.pRMS(2,:) = sqrt(est.RMS(2,:)/(N-1));                 %KF Root Mean Square Estimation Error of x2

est.Infmean(1,:) = est.Infmean(1,:)/N;                      %IF Estimation Error Sample Mean of x1
est.Infmean(2,:) = est.Infmean(2,:)/N;                      %IF Estimation Error Sample Mean of x2

est.pInfRMS(1,:) = sqrt(est.InfRMS(1,:)/(N-1));                 %IF Root Mean Square Estimation Error of x1
est.pInfRMS(2,:) = sqrt(est.InfRMS(2,:)/(N-1));                 %IF Root Mean Square Estimation Error of x2

figure(2)
plot(t,est.pRMS(1,:),'k','LineStyle','--','LineWidth',2);
plot(t,-est.pRMS(1,:),'k','LineStyle','--','LineWidth',2);
plot(t,est.mean(1,:),'b','LineStyle','-','LineWidth',2)
plot(t,est.Pplot(1,:),'r','LineStyle','-','LineWidth',2)
plot(t,-est.Pplot(1,:),'r','LineStyle','-','LineWidth',2)

figure(3)
plot(t,est.pRMS(2,:),'k','LineStyle','--','LineWidth',2);
plot(t,-est.pRMS(2,:),'k','LineStyle','--','LineWidth',2);
plot(t,est.mean(2,:),'b','LineStyle','-','LineWidth',2)
plot(t,est.Pplot(2,:),'r','LineStyle','-','LineWidth',2)
plot(t,-est.Pplot(2,:),'r','LineStyle','-','LineWidth',2)

figure(4)
plot(t,est.pInfRMS(1,:),'k','LineStyle','--','LineWidth',2);
plot(t,-est.pInfRMS(1,:),'k','LineStyle','--','LineWidth',2);
plot(t,est.Infmean(1,:),'b','LineStyle','-','LineWidth',2)
plot(t,est.InfPplot(1,:),'r','LineStyle','-','LineWidth',2)
plot(t,-est.InfPplot(1,:),'r','LineStyle','-','LineWidth',2)

figure(5)
plot(t,est.pInfRMS(2,:),'k','LineStyle','--','LineWidth',2);
plot(t,-est.pInfRMS(2,:),'k','LineStyle','--','LineWidth',2);
plot(t,est.Infmean(2,:),'b','LineStyle','-','LineWidth',2)
plot(t,est.InfPplot(2,:),'r','LineStyle','-','LineWidth',2)
plot(t,-est.InfPplot(2,:),'r','LineStyle','-','LineWidth',2)

figure(1)
plot(t,y_bar,'--');         %Adding y_bar to figure 1