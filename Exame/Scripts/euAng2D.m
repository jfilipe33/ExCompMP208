function Dout = euAng2D(a0)
    phi = a0(1);
    theta = a0(2);
    psi  = a0(3);
    
    Dout = [cos(psi)*cos(theta) (cos(psi)*sin(theta)*sin(phi))+(sin(psi)*cos(phi)) -(cos(psi)*sin(theta)*cos(phi))+(sin(psi)*sin(phi));
        -sin(psi)*cos(theta) -(sin(psi)*sin(theta)*sin(phi))+(cos(psi)*cos(phi)) (sin(psi)*sin(theta)*cos(phi))+(cos(psi)*sin(phi));
        sin(theta) -cos(theta)*sin(phi) cos(theta)*cos(phi)];
    
    %obj.Omegabar = [cos(theta)*cos(psi) sin(psi) 0; -cos(theta)*sin(psi) cos(psi) 0; sin(theta) 0 1];
        
end
