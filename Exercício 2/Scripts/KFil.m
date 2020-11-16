%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%              - Exercício Computacional de MP208 -             %
%    --- Optimal Filtering with Aerospace Applications ---      %
%                                                               %
%              Autor: João Filipe R. P de A. Silva              %
%                                                               %
%               Function Script: Filtro de Kalman               %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function est = KFil(sys,est)

%Prediction
    x_k1k = sys.A*est.x + sys.B*sys.u;      %Predictive States Estimative
    y_k1k = sys.C*x_k1k;                    %Predictive Measure

    P_k1k = sys.A*est.P*sys.A' + est.Q;     %Predictive States Estimative Covariance
    PY_k1k = sys.C*P_k1k*sys.C'+ est.R;     %Predictive Measure Covariance
    PXY_k1k = P_k1k*sys.C';                 %Predictive Measure Cross-Covariance
    
%Update
    K_k1 = PXY_k1k*inv(PY_k1k);                     %Kalman Filter Gain
    x_k1k1 = x_k1k + K_k1*(sys.y - y_k1k);          %Filtered Estimative
    P_k1k1 = P_k1k - PXY_k1k*inv(PY_k1k)*PXY_k1k';  %Estimation Error Covariance
    
%Updated Variables
    est.x = x_k1k1;
    est.P = P_k1k1;
    
    end