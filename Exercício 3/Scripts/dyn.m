function xnew = dyn(x,sys)         
    
    xnew = zeros(2,1);
    xnew(1) = -x(1) + x(2);
    xnew(2) = -0.1*x(1)^2 - 1 + sys.u;
    xnew = xnew + sys.w;
    
end