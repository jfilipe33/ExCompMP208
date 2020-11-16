%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%              - Exercício Computacional de MP208 -             %
%    --- Optimal Filtering with Aerospace Applications ---      %
%                                                               %
%              Autor: João Filipe R. P de A. Silva              %
%                                                               %
%             Function Script: Filtro de Informação             %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function est = InfFil(sys,est)

    est.L = inv(est.P);                 %Filtered Information Matrix
    est.z = est.L*est.x;                %Information Filter Estimated States

%Prediction

    Pi_k = inv(sys.A')*est.L*inv(sys.A);        %Auxiliary Calculation Matrix
    K_k = Pi_k*inv(Pi_k + inv(est.Q));          %Auxiliary Calculation Matrix 
        
    z_k1k = (eye(2) - K_k)*inv(sys.A)'*est.z + (eye(2) - K_k)*Pi_k*sys.B*sys.u; %Predictive Information States Estimative
    L_k1k = (eye(2) - K_k)*Pi_k;                                                %Predictive States Estimative Information
    
%Update
    
    z_k1k1 = z_k1k + sys.C'*inv(est.R)*sys.y;   %Filtered Estimative
    L_k1k1 = L_k1k + sys.C'*inv(est.R)*sys.C;   %Estimation Error Information
    
%Updated Variables
    est.z = z_k1k1;
    est.L = L_k1k1;
    
    end