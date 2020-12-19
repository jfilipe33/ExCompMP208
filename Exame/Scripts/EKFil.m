%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%              - Exercício Computacional de MP208 -             %
%    --- Optimal Filtering with Aerospace Applications ---      %
%                                                               %
%              Autor: João Filipe R. P de A. Silva              %
%                                                               %
%          Function Script: Filtro de Kalman Estendido          %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ekf = EKFil(sys,ekf,simul)

%Prediction
    
       
    estimates = ekf.x;
    auxP = ekf.P;
    
    %integration using RK4
    
    ekf.F = matF(estimates,sys);    %Jacobian of function f
    ekf.G = matG(estimates);        
    k1x = simul.Ts*dynFil(estimates,sys);
    k1P = simul.Ts*dynCov(auxP,sys,ekf);
    
    ekf.F = matF(estimates+k1x/2,sys);  
    ekf.G = matG(estimates+k1x/2);
    k2x = simul.Ts*dynFil(estimates+k1x/2,sys);
    k2P = simul.Ts*dynCov(auxP+k1P/2,sys,ekf);
    
    ekf.F = matF(estimates+k2x/2,sys);  
    ekf.G = matG(estimates+k2x/2);
    k3x = simul.Ts*dynFil(estimates+k2x/2,sys); 
    k3P = simul.Ts*dynCov(auxP+k2P/2,sys,ekf);
    
    ekf.F = matF(estimates+k3x,sys);  
    ekf.G = matG(estimates+k3x);
    k4x = simul.Ts*dynFil(estimates+k3x,sys); 
    k4P = simul.Ts*dynCov(auxP+k3P,sys,ekf);
    estimates  = estimates + k1x/6 + k2x/3 + k3x/3 + k4x/6; 
    auxP = auxP + k1P/6 + k2P/3 + k3P/3 + k4P/6;
    
    x_k1k = estimates;         %Predictive States Estimative
    
    ekf.H = matH(x_k1k,sys);    %Jacobian of function h
    
    ekf.D = euAng2D(x_k1k(1:3));
    ekf.h = [ekf.D*sys.e3; ekf.D*sys.n];
    y_k1k = ekf.h;             %Predictive Measure

    P_k1k = auxP;                            %Predictive States Estimative Covariance
    PY_k1k = ekf.H*P_k1k*ekf.H' + ekf.R;     %Predictive Measure Covariance
    PXY_k1k = P_k1k*ekf.H';                  %Predictive Measure Cross-Covariance
    
%Update
    
    K_k1 = PXY_k1k*inv(PY_k1k);                     %Extended Kalman Filter Gain
    x_k1k1 = x_k1k + K_k1*(sys.y - y_k1k);          %Filtered Estimative
    P_k1k1 = P_k1k - PXY_k1k*inv(PY_k1k)*PXY_k1k';  %Estimation Error Covariance
    
%Updated Variables
    ekf.x = x_k1k1;
    ekf.P = P_k1k1;
    
    end