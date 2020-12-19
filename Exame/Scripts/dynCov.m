function Pnew = dynCov(P,sys,ekf)         
    
    Pnew = zeros(6,6);
    Pnew = ekf.F*P + P*ekf.F' + ekf.G*ekf.Q*ekf.G';
    
end