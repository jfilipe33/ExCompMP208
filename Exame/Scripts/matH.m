function H = matH(x,sys)
       
    phi = x(1);
    theta = x(2);
    psi = x(3);

    dDphi = [0 (cos(psi)*sin(theta)*cos(phi) - sin(psi)*sin(phi)) (cos(psi)*sin(theta)*sin(phi) + sin(psi)*cos(phi));
             0 (-sin(psi)*sin(theta)*cos(phi) - cos(psi)*sin(phi)) (-sin(psi)*sin(theta)*sin(phi) + cos(psi)*cos(phi));
             0 -cos(theta)*cos(phi) -cos(theta)*sin(phi)];
    dDtheta = [-cos(psi)*sin(theta) cos(psi)*cos(theta)*sin(phi) -cos(psi)*cos(theta)*cos(phi);
               sin(psi)*sin(theta) -sin(psi)*cos(theta)*sin(phi) sin(psi)*cos(theta)*cos(phi);
               cos(theta) sin(theta)*sin(phi) -sin(theta)*cos(phi)];
    dDpsi = [-sin(psi)*cos(theta) (-sin(psi)*sin(theta)*sin(phi) + cos(psi)*cos(phi)) (sin(psi)*sin(theta)*cos(phi) + cos(psi)*sin(phi));
             -cos(psi)*cos(theta) (-cos(psi)*sin(theta)*sin(phi) - sin(psi)*cos(phi)) (cos(psi)*sin(theta)*cos(phi) - sin(psi)*sin(phi));
             0 0 0];
    
    H = [dDphi*sys.e3 dDtheta*sys.e3 dDpsi*sys.e3 zeros(3,3);
         dDphi*sys.n dDtheta*sys.n dDpsi*sys.n zeros(3,3)];
end