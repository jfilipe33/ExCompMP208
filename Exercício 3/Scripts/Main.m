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

sys.R = 0.01;                        %Measurement Noise variance
sys.Q = diag([1E-2 1E-2]);           %States Disturbance variance
sys.P = diag([1E0 1E0]);            %X Covariance Matrix
sys.m = [0;0];                      %State Expected Value
sys.Ts = 0.1;                       %Sampling Time
t = 0:sys.Ts:50;                   %Simulation Period
N = 100;                            %Number of Realizations in the simulation
est.mean = zeros(2,length(t));          %Error Sample Mean 
sys.mean = zeros(1,length(t));          %Error Sample Mean 
est.RMS = zeros(2,length(t));           %RMS Error 
est.sig = zeros(2,length(t));           %Estimation Error SD Mean
est.Pkk = zeros(2,2);


%% Simulation
  
for k = 1:N
    sys.x = sqrtm(sys.P)*randn(2,1) + sys.m;             %Initial value of X
    sys.y = sys.x(1);                                    %Output Initial Value
    est.P = diag([1E0 1E0]);                             %Initial Estimate Covariance Matrix
    est.x = sys.m;                                       %Initial Estimate of X
    est.Q = diag([1E-5 5E-3]);
    est.R = 0.01;

    for cont = 1:length(t)
    
        %Model Behaviour
        
        sys.u = -10*sys.y + 10; 
        states = sys.x;
        sys.w = sqrtm(sys.Q)*randn(2,1);                %State Disturbance Realization
        sys.v = sqrt(sys.R)*randn(1,1);                 %Measurement Noise Realization

        %integration using RK4

        k1 = sys.Ts*dyn(states,sys);
        k2 = sys.Ts*dyn(states+k1/2,sys);
        k3 = sys.Ts*dyn(states+k2/2,sys);         
        k4 = sys.Ts*dyn(states+k3,sys); 
        states  = states + k1/6 + k2/3 + k3/3 + k4/6; 

        sys.x = states;
        sys.y = sys.x(1) + sys.v;

        % Extended Kalman Filter   
        
        est = EKFil(sys,est);
        ex = est.x - sys.x;                             %Extended Kalman Filter estimation Error
        
        ex1(k,cont) = ex(1);                            %Isolating ex component 1
        ex2(k,cont) = ex(2);                            %Isolating ex component 2
        
        yp(k,cont) = sys.y;                                 %Output plotting variable;
        
        sys.mean(cont) = sys.mean(cont) + sys.y;            %Sample mean of y
        
        est.mean(1,cont) = est.mean(1,cont) + ex(1);          %Sample mean of EKF Estimation Error of x1
        est.mean(2,cont) = est.mean(1,cont) + ex(2);          %Sample mean of EKF Estimation Error of x2

        est.RMS(1,cont) = est.RMS(1,cont) + ex(1)*ex(1)';     %Sum of Squared Errors in ex(1)
        est.RMS(2,cont) = est.RMS(2,cont) + ex(2)*ex(2)';     %Sum of Squared Errors in ex(2)
        
        est.Pkk = sqrtm(est.P);
        est.sig(1,cont) = est.Pkk(1,1);          %Theoretical EKF Standard Deviation of ex(1)
        est.sig(2,cont) = est.Pkk(2,2);          %Theoretical EKF Standard Deviation of ex(1)

    end
    
    figure(1)
    plot(t,ex1(k,:));
    hold on;
    title('\textbf{Extended Kalman Filter Estimation Error of state $x_1$}','interpreter','LaTeX');
    xlabel('\textbf{Time [s]}','interpreter','LaTeX');
    ylabel('\textbf{Amplitude [\,-\,]}','interpreter','LaTeX');
    
    figure(2)
    plot(t,ex2(k,:));
    hold on;
    title('\textbf{Extended Kalman Filter Estimation Error of state $x_2$}','interpreter','LaTeX');
    xlabel('\textbf{Time [s]}','interpreter','LaTeX');
    ylabel('\textbf{Amplitude [\,-\,]}','interpreter','LaTeX');
    
    figure(3)
    plot(t,yp(k,:));
    hold on;
    title('\textbf{Measured Outputs $y_k$}','interpreter','LaTeX');
    xlabel('\textbf{Time [s]}','interpreter','LaTeX');
    ylabel('\textbf{Amplitude [\,-\,]}','interpreter','LaTeX');
    
end

est.mean(1,:) = est.mean(1,:)/N;                      %EKF Estimation Error Sample Mean of x1
est.mean(2,:) = est.mean(2,:)/N;                      %EKF Estimation Error Sample Mean of x2

sys.mean = sys.mean/N;                                %EKF Sample Mean of y

est.pRMS(1,:) = sqrt(est.RMS(1,:)/N);                 %EKF Root Mean Square Estimation Error of x1
est.pRMS(2,:) = sqrt(est.RMS(2,:)/N);                 %EKF Root Mean Square Estimation Error of x2

figure(1)
plot(t,est.pRMS(1,:),'k','LineStyle','--','LineWidth',2);
plot(t,-est.pRMS(1,:),'k','LineStyle','--','LineWidth',2);
plot(t,est.mean(1,:),'y','LineStyle','-','LineWidth',2)
plot(t,est.sig(1,:),'r','LineStyle','-','LineWidth',2)
plot(t,-est.sig(1,:),'r','LineStyle','-','LineWidth',2)

figure(2)
plot(t,est.pRMS(2,:),'k','LineStyle','--','LineWidth',2);
plot(t,-est.pRMS(2,:),'k','LineStyle','--','LineWidth',2);
plot(t,est.mean(2,:),'y','LineStyle','-','LineWidth',2)
plot(t,est.sig(2,:),'r','LineStyle','-','LineWidth',2)
plot(t,-est.sig(2,:),'r','LineStyle','-','LineWidth',2)

figure(3)
plot(t,sys.mean,'y','LineStyle','-','LineWidth',2)
