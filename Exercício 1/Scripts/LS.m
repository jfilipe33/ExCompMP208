%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%              - Exercício Computacional de MP208 -             %
%    --- Optimal Filtering with Aerospace Applications ---      %
%                                                               %
%              Autor: João Filipe R. P de A. Silva              %
%                                                               %
%           Function Script: Estimador Least Squares            %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ls = LS(sys,ls,n)
    
    ls.A = zeros(2,2);              %Auxiliary LS Sum Matrix
    ls.B = zeros(2,1);              %Auxiliary LS Sum Matrix
    
    for m = 1:sys.N                                     %Measurement Batch Loop
        v = sys.R*randn(1);                             %Measurement Noise 
        sys.y(m,1) = sys.H(m,:)*sys.theta + v;          %Measurement Equation
        ls.A = ls.A + sys.H(m,:)'*sys.W*sys.H(m,:);
        ls.B = ls.B + sys.H(m,:)'*sys.W*sys.y(m,1);
    end
    
    ls.thetaLS(:,n) = inv(ls.A)*(ls.B);
    ls.RMS_LS = ls.RMS_LS + ((ls.thetaLS(:,n)-sys.theta)*(ls.thetaLS(:,n)-sys.theta)');

end