%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%              - Exercício Computacional de MP208 -             %
%    --- Optimal Filtering with Aerospace Applications ---      %
%                                                               %
%              Autor: João Filipe R. P de A. Silva              %
%                                                               %
%          Function Script: Filtro de Kalman Unscented          %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ukf = UKFil(sys,ukf,simul)

    ukf.xaug = [ukf.x; zeros(6,1); zeros(6,1)];
    ukf.Paug = blkdiag(ukf.P,ukf.Q,ukf.R);
    
    x_k1k = zeros(6,1);
    y_k1k = zeros(6,1);
    P_k1k = zeros(6,6);
    PY_k1k = zeros(6,6);
    PXY_k1k = zeros(6,6);
    
    ukf = SP(ukf);              %Sigma Points calculation
      
    for j = 1:length(ukf.SP')
        estimates = ukf.SP((1:6),j);
        W = ukf.SP((7:12),j);
        V = ukf.SP((13:18),j);
        k1 = simul.Ts*dynUFil(estimates,W,sys);
        k2 = simul.Ts*dynUFil(estimates+k1/2,W,sys);
        k3 = simul.Ts*dynUFil(estimates+k2/2,W,sys);         
        k4 = simul.Ts*dynUFil(estimates+k3,W,sys); 
        X_k1k(:,j) = estimates + k1/6 + k2/3 + k3/3 + k4/6;
        
        ukf.D = euAng2D(X_k1k((1:3),j));
        ukf.h = [ukf.D*sys.e3; ukf.D*sys.n];
        Y_k1k(:,j) = ukf.h + V;                
        
        x_k1k = x_k1k + ukf.ro(j)*X_k1k(:,j);   %Predictive States 
        y_k1k = y_k1k + ukf.ro(j)*Y_k1k(:,j);   %Predictive Measure
    end
    
%Prediction
     
    for j = 1:length(ukf.SP')
        P_k1k = P_k1k + ukf.ro(j)*(X_k1k(:,j) - x_k1k)*(X_k1k(:,j) - x_k1k)';          %Predictive States Estimative Covariance
        PY_k1k = PY_k1k + ukf.ro(j)*(Y_k1k(:,j) - y_k1k)*(Y_k1k(:,j) - y_k1k)';        %Predictive Measure Covariance
        PXY_k1k = PXY_k1k + ukf.ro(j)*(X_k1k(:,j) - x_k1k)*(Y_k1k(:,j) - y_k1k)';      %Predictive Measure Cross-Covariance
    end

%Update
    
    K_k1 = PXY_k1k*inv(PY_k1k);                     %Extended Kalman Filter Gain
    x_k1k1 = x_k1k + K_k1*(sys.y - y_k1k);          %Filtered Estimative
    P_k1k1 = P_k1k - PXY_k1k*inv(PY_k1k)*PXY_k1k';  %Estimation Error Covariance
    
%Updated Variables
    ukf.x = x_k1k1(1:6);
    ukf.P = P_k1k1(1:6,1:6);
    
end