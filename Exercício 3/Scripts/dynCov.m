function Pnew = dynCov(P,sys,est)         
    
    Pnew = zeros(2,2);
    Pnew = est.F*P + P*est.F' + est.Q;
    
end