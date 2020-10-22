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

function ls = LS(sys,ls,i)
    
    ls.A = ls.A + sys.H'*sys.W*sys.H;
    ls.B = ls.B + sys.H'*sys.W*sys.y;
    ls.thetaLS(:,i) = inv(ls.A)*(ls.B);
    ls.RMS_LS = ls.RMS_LS + ((ls.thetaLS(:,i)-sys.theta)*(ls.thetaLS(:,i)-sys.theta)');

end