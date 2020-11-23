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

function est = EKFil(sys,est)

%Prediction
    est.F = [-1 1; -0.2*est.x(1) 0];
    est.H = [1 0];
    estimates = est.x;
    auxP = est.P;
    
    %integration using RK4
    
    k1s = sys.Ts*dyn(estimates,sys);
    k1c = sys.Ts*dynCov(auxP,sys,est);
    k2s = sys.Ts*dyn(estimates+k1s/2,sys);
    k2c = sys.Ts*dynCov(auxP+k1c/2,sys,est);
    k3s = sys.Ts*dyn(estimates+k2s/2,sys); 
    k3c = sys.Ts*dynCov(auxP+k2c/2,sys,est);
    k4s = sys.Ts*dyn(estimates+k3s,sys); 
    k4c = sys.Ts*dynCov(auxP+k3c,sys,est);
    estimates  = estimates + k1s/6 + k2s/3 + k3s/3 + k4s/6; 
    auxP = auxP + k1c/6 + k2c/3 + k3c/3 + k4c/6;
    
    x_k1k = estimates;         %Predictive States Estimative
    y_k1k = estimates(1);      %Predictive Measure

    P_k1k = auxP;                            %Predictive States Estimative Covariance
    PY_k1k = est.H*P_k1k*est.H' + est.R;     %Predictive Measure Covariance
    PXY_k1k = P_k1k*est.H';                  %Predictive Measure Cross-Covariance
    
%Update
    K_k1 = PXY_k1k*inv(PY_k1k);                     %Extended Kalman Filter Gain
    x_k1k1 = x_k1k + K_k1*(sys.y - y_k1k);          %Filtered Estimative
    P_k1k1 = P_k1k - PXY_k1k*inv(PY_k1k)*PXY_k1k';  %Estimation Error Covariance
    
%Updated Variables
    est.x = x_k1k1;
    est.P = P_k1k1;
    
    end