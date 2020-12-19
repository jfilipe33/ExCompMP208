function xnew = dyn(x,sys)         
    
    xnew = zeros(6,1);
    euAng = x(1:3); 
    
    alphadot = A(euAng)*sys.omega;
    betadot = sys.w(4:6);
    xnew = [alphadot;betadot];
    
end