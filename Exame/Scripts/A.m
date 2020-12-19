function mat = A(alpha)         
    
    mat = [cos(alpha(3))/cos(alpha(2)) -(sin(alpha(3))/cos(alpha(2))) 0;
           sin(alpha(3)) cos(alpha(3)) 0;
           -(cos(alpha(3))*sin(alpha(2))/cos(alpha(2))) sin(alpha(3))*sin(alpha(2))/cos(alpha(2)) 1];
    
end