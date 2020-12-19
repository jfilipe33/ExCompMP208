function G = matG(x)

    euAng = x(1:3);
    G = [-A(euAng) zeros(3,3);
        zeros(3,3) eye(3)];

end