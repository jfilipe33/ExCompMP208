function F = matF(x,sys)
       
    u = sys.u - x(4:6);
    phi = x(1);
    theta = x(2);
    psi = x(3);

    dAphi = zeros(3,3);
    dAtheta = [((cos(psi)*sin(theta))/(cos(theta)^2)) ((-sin(psi)*sin(theta))/(cos(theta)^2)) 0;
               0 0 0;
               -cos(psi)/(cos(theta)^2) sin(psi)/(cos(theta)^2) 0];
    dApsi = [-sin(psi)/cos(theta) -cos(psi)/cos(theta) 0;
             cos(psi) -sin(psi) 0;
             (sin(psi)*sin(theta))/cos(theta) (cos(psi)*sin(theta))/cos(theta) 0];
    
    F = [dAphi*u dAtheta*u dApsi*u -A([phi theta psi]);
         zeros(3,6)];
end