function xnew = dynFil(x,sys)         
    
    
    xnew = zeros(6,1);
    euAng = x(1:3);
    bias = x(4:6); 
    u = sys.u - bias;
    
    xnew = [A(euAng)*u;
            zeros(3,1)];
    
end