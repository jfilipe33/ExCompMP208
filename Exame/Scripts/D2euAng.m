function ang = D2euAng(D)
      
    phi = -atan2(D(3,2),D(3,3));
    theta = asin(D(3,1));
    psi = -atan2(D(2,1),D(1,1));
    ang = [phi theta psi]';
        
end
