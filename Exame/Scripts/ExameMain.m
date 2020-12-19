%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%              - Exercício Computacional de MP208 -             %
%    --- Optimal Filtering with Aerospace Applications ---      %
%                                                               %
%              Autor: João Filipe R. P de A. Silva              %
%                                                               %
%      Main Script: Exame Final - Determinação de Atitude       %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cache Clean-up
clear all
close all
clc

%% System parameters

tstart = clock;                     %Time counter initialization
sys.Qgy = 2.5E-7 * eye(3);          %Rate-Gyro Noise Covariance
sys.Qdgy = 1E-12 * eye(3);          %Rate-Gyro Drift Covariance
sys.Qac = 1.5E-5 * eye(3);          %Accelerometer Covariance
sys.Qmg = 0.02 * eye(3);            %Magnetometer Covariance
sys.mg = [13.7 -4.6 -10.9]';        %Local Magnetic Field
sys.g = 9.81;                       %Gravity Acceleration
sys.P = diag([((pi/9)^2)*ones(1,3) 0.1*ones(1,3)]);     %States Covariance Matrix
sys.m = zeros(6,1);                 %States expected value
sys.e3 = [0 0 1]';                  %Unit Vector
sys.n = sys.mg/norm(sys.mg);        %Magnectic Field unit Vector

%% Simulation Parameters

simul.Ts = 0.01;
simul.t = 0:simul.Ts:10-simul.Ts;

ekf.sig = zeros(6,length(simul.t));           %Estimation Error SD Mean
ekf.Pkk = zeros(6,6);                           

ukf.sig = zeros(6,length(simul.t));           %Estimation Error SD Mean
ukf.Pkk = zeros(6,6);


%% Simulation


for N = 1:100
    
    sys.x = sqrtm(sys.P)*randn(6,1) + sys.m;        %States initialization
    
    % EKF parameters

    ekf.Q = blkdiag(1E0*sys.Qgy,sys.Qdgy);
    ekf.R = blkdiag(1E2*sys.Qac,1E-1*sys.Qmg);
    ekf.P = 1E0*sys.P;
    ekf.x = sys.m;
    
    % UKF Parameters
    
    ukf.Q = blkdiag(1E0*sys.Qgy,sys.Qdgy);
    ukf.R = blkdiag(1E2*sys.Qac,1E-1*sys.Qmg);
    ukf.P = sys.P;
    ukf.x = sys.m;
    

    for k = 1:length(simul.t)

        %Angular Velocity Reference
        sys.omega = (2/20) * pi * [sin(0.5*pi*simul.Ts*k) sin(0.5*pi*simul.Ts*k + pi/2) sin(0.5*pi*simul.Ts*k + pi)]';

        %Model Behaviour

        states = sys.x;
        sys.w = [sqrtm(sys.Qgy)*randn(3,1); sqrtm(sys.Qdgy)*randn(3,1)];                         %State Disturbance Realization
        sys.v = [sqrtm(sys.Qac)*randn(3,1)/sys.g; sqrtm(sys.Qmg)*randn(3,1)/norm(sys.mg)];       %Measurement Noise Realization
        sys.omegac = sys.omega + states(4:6) + sys.w(1:3);
        sys.u = sys.omegac;

        %integration using RK4

        k1 = simul.Ts*dyn(states,sys);
        k2 = simul.Ts*dyn(states+k1/2,sys);
        k3 = simul.Ts*dyn(states+k2/2,sys);         
        k4 = simul.Ts*dyn(states+k3,sys); 
        states = states + k1/6 + k2/3 + k3/3 + k4/6;

        %Updating states

        sys.x = states;
        sys.alpha = states(1:3);
        sys.beta = states(4:6);

        %Sensor measurements

        sys.D = euAng2D(sys.alpha);
        sys.h = [sys.D*sys.e3;sys.D*sys.n];
        sys.y = sys.h + sys.v;
        
        % Extended Kalman Filter   
        
%         ekf = EKFil(sys,ekf,simul);
%         ex = ekf.x - sys.x;                             %Extended Kalman Filter estimation Error
%         
%         ex1(N,k) = ex(1);                            %Isolating ex component 1
%         ex2(N,k) = ex(2);                            %Isolating ex component 2
%         ex3(N,k) = ex(3);                            %Isolating ex component 3
%         ex4(N,k) = ex(4);                            %Isolating ex component 4
%         ex5(N,k) = ex(5);                            %Isolating ex component 5
%         ex6(N,k) = ex(6);                            %Isolating ex component 6
%         
%         
%         ekf.Pkk = sqrtm(ekf.P);
%         ekf.sig(1,k) = ekf.Pkk(1,1);          %Theoretical EKF Standard Deviation of ex(1)
%         ekf.sig(2,k) = ekf.Pkk(2,2);          %Theoretical EKF Standard Deviation of ex(2)
%         ekf.sig(3,k) = ekf.Pkk(3,3);          %Theoretical EKF Standard Deviation of ex(3)
%         ekf.sig(4,k) = ekf.Pkk(4,4);          %Theoretical EKF Standard Deviation of ex(4)
%         ekf.sig(5,k) = ekf.Pkk(5,5);          %Theoretical EKF Standard Deviation of ex(5)
%         ekf.sig(6,k) = ekf.Pkk(6,6);          %Theoretical EKF Standard Deviation of ex(6)
        
        % Unscented Kalman Filter 
        
        ukf = UKFil(sys,ukf,simul);
        ux = ukf.x - sys.x;
        
        ux1(N,k) = ux(1);                            %Isolating ux component 1
        ux2(N,k) = ux(2);                            %Isolating ux component 2
        ux3(N,k) = ux(3);                            %Isolating ux component 3
        ux4(N,k) = ux(4);                            %Isolating ux component 4
        ux5(N,k) = ux(5);                            %Isolating ux component 5
        ux6(N,k) = ux(6);                            %Isolating ux component 6
        
        ukf.Pkk = sqrtm(ukf.P);
        ukf.sig(1,k) = ukf.Pkk(1,1);          %Theoretical UKF Standard Deviation of ux(1)
        ukf.sig(2,k) = ukf.Pkk(2,2);          %Theoretical UKF Standard Deviation of ux(2)
        ukf.sig(3,k) = ukf.Pkk(3,3);          %Theoretical UKF Standard Deviation of ux(3)
        ukf.sig(4,k) = ukf.Pkk(4,4);          %Theoretical UKF Standard Deviation of ux(4)
        ukf.sig(5,k) = ukf.Pkk(5,5);          %Theoretical UKF Standard Deviation of ux(5)
        ukf.sig(6,k) = ukf.Pkk(6,6);          %Theoretical UKF Standard Deviation of ux(6)
            
        %Plotting history

        histOmega(:,k) = sys.omega;
        histOmegaC(:,k) = sys.omegac;

    end
    
    
%     figure(1)
%     plot(simul.t,rad2deg(ex1(N,:)))
%     hold on;
%     title('\textbf{Extended Kalman Filter Estimation Error of state $x_1$}','interpreter','LaTeX');
%     xlabel('\textbf{Time [s]}','interpreter','LaTeX');
%     ylabel('\textbf{Angle [deg]}','interpreter','LaTeX');
% 
%     figure(2)
%     plot(simul.t,rad2deg(ex2(N,:)))
%     hold on;
%     title('\textbf{Extended Kalman Filter Estimation Error of state $x_2$}','interpreter','LaTeX');
%     xlabel('\textbf{Time [s]}','interpreter','LaTeX');
%     ylabel('\textbf{Angle [deg]}','interpreter','LaTeX');
% 
%     figure(3)
%     plot(simul.t,rad2deg(ex3(N,:)))
%     hold on;
%     title('\textbf{Extended Kalman Filter Estimation Error of state $x_3$}','interpreter','LaTeX');
%     xlabel('\textbf{Time [s]}','interpreter','LaTeX');
%     ylabel('\textbf{Angle [deg]}','interpreter','LaTeX');
%     
%     figure(4)
%     plot(simul.t,ex4(N,:))
%     hold on
%     title('\textbf{Extended Kalman Filter Estimation Error of state $x_4$}','interpreter','LaTeX');
%     xlabel('\textbf{Time [s]}','interpreter','LaTeX');
%     ylabel('\textbf{Angle [deg]}','interpreter','LaTeX');
% 
%     figure(5)
%     plot(simul.t,ex5(N,:))
%     hold on
%     title('\textbf{Extended Kalman Filter Estimation Error of state $x_5$}','interpreter','LaTeX');
%     xlabel('\textbf{Time [s]}','interpreter','LaTeX');
%     ylabel('\textbf{Angle [deg]}','interpreter','LaTeX');
% 
%     figure(6)
%     plot(simul.t,ex6(N,:))
%     hold on
%     title('\textbf{Extended Kalman Filter Estimation Error of state $x_6$}','interpreter','LaTeX');
%     xlabel('\textbf{Time [s]}','interpreter','LaTeX');
%     ylabel('\textbf{Angle [deg]}','interpreter','LaTeX');
    
    figure(1)
    plot(simul.t,rad2deg(ux1(N,:)))
    hold on;
    title('\textbf{Unscented Kalman Filter Estimation Error of state $x_1$}','interpreter','LaTeX');
    xlabel('\textbf{Time [s]}','interpreter','LaTeX');
    ylabel('\textbf{Angle [deg]}','interpreter','LaTeX');

    figure(2)
    plot(simul.t,rad2deg(ux2(N,:)))
    hold on;
    title('\textbf{Unscented Kalman Filter Estimation Error of state $x_2$}','interpreter','LaTeX');
    xlabel('\textbf{Time [s]}','interpreter','LaTeX');
    ylabel('\textbf{Angle [deg]}','interpreter','LaTeX');

    figure(3)
    plot(simul.t,rad2deg(ux3(N,:)))
    hold on;
    title('\textbf{Unscented Kalman Filter Estimation Error of state $x_3$}','interpreter','LaTeX');
    xlabel('\textbf{Time [s]}','interpreter','LaTeX');
    ylabel('\textbf{Angle [deg]}','interpreter','LaTeX');
    
    figure(4)
    plot(simul.t,ux4(N,:))
    hold on
    title('\textbf{Unscented Kalman Filter Estimation Error of state $x_4$}','interpreter','LaTeX');
    xlabel('\textbf{Time [s]}','interpreter','LaTeX');
    ylabel('\textbf{Angle [deg]}','interpreter','LaTeX');

    figure(5)
    plot(simul.t,ux5(N,:))
    hold on
    title('\textbf{Unscented Kalman Filter Estimation Error of state $x_5$}','interpreter','LaTeX');
    xlabel('\textbf{Time [s]}','interpreter','LaTeX');
    ylabel('\textbf{Angle [deg]}','interpreter','LaTeX');

    figure(6)
    plot(simul.t,ux6(N,:))
    hold on
    title('\textbf{Unscented Kalman Filter Estimation Error of state $x_6$}','interpreter','LaTeX');
    xlabel('\textbf{Time [s]}','interpreter','LaTeX');
    ylabel('\textbf{Angle [deg]}','interpreter','LaTeX');
end

% figure(1)
% plot(simul.t,ekf.sig(1,:),'r','LineStyle','--','LineWidth',2.5)
% plot(simul.t,-ekf.sig(1,:),'r','LineStyle','--','LineWidth',2.5)
% figure(2)
% plot(simul.t,ekf.sig(2,:),'r','LineStyle','--','LineWidth',2.5)
% plot(simul.t,-ekf.sig(2,:),'r','LineStyle','--','LineWidth',2.5)
% figure(3)
% plot(simul.t,ekf.sig(3,:),'r','LineStyle','--','LineWidth',2.5)
% plot(simul.t,-ekf.sig(3,:),'r','LineStyle','--','LineWidth',2.5)
% figure(4)
% plot(simul.t,ekf.sig(4,:),'r','LineStyle','--','LineWidth',2.5)
% plot(simul.t,-ekf.sig(4,:),'r','LineStyle','--','LineWidth',2.5)
% figure(5)
% plot(simul.t,ekf.sig(5,:),'r','LineStyle','--','LineWidth',2.5)
% plot(simul.t,-ekf.sig(5,:),'r','LineStyle','--','LineWidth',2.5)
% figure(6)
% plot(simul.t,ekf.sig(6,:),'r','LineStyle','--','LineWidth',2.5)
% plot(simul.t,-ekf.sig(6,:),'r','LineStyle','--','LineWidth',2.5)

figure(1)
plot(simul.t,ukf.sig(1,:),'g','LineStyle','--','LineWidth',2.5)
plot(simul.t,-ukf.sig(1,:),'g','LineStyle','--','LineWidth',2.5)
figure(2)
plot(simul.t,ukf.sig(2,:),'g','LineStyle','--','LineWidth',2.5)
plot(simul.t,-ukf.sig(2,:),'g','LineStyle','--','LineWidth',2.5)
figure(3)
plot(simul.t,ukf.sig(3,:),'g','LineStyle','--','LineWidth',2.5)
plot(simul.t,-ukf.sig(3,:),'g','LineStyle','--','LineWidth',2.5)
figure(4)
plot(simul.t,ukf.sig(4,:),'g','LineStyle','--','LineWidth',2.5)
plot(simul.t,-ukf.sig(4,:),'g','LineStyle','--','LineWidth',2.5)
figure(5)
plot(simul.t,ukf.sig(5,:),'g','LineStyle','--','LineWidth',2.5)
plot(simul.t,-ukf.sig(5,:),'g','LineStyle','--','LineWidth',2.5)
figure(6)
plot(simul.t,ukf.sig(6,:),'g','LineStyle','--','LineWidth',2.5)
plot(simul.t,-ukf.sig(6,:),'g','LineStyle','--','LineWidth',2.5)


figure(7)
plot(simul.t,histOmega(1,:),'--m',simul.t,histOmega(2,:),'--k',simul.t,histOmega(3,:),'--g')
hold on
plot(simul.t,histOmegaC(1,:),'m',simul.t,histOmegaC(2,:),'k',simul.t,histOmegaC(3,:),'g')
title('\textbf{True Omega $\omega^{b/g}_b$ vs. Measured Omega $\check{\omega}_b$}','interpreter','LaTeX');
xlabel('\textbf{Time [s]}','interpreter','LaTeX');
ylabel('\textbf{Angular Velocity [deg/s]}','interpreter','LaTeX');

etime(clock,tstart) %Calculate simulation elapsed time
