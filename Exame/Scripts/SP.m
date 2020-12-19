function ukf = SP(ukf)

    n = length(ukf.xaug);
    kappa = 3-n;
    S = sqrtm(ukf.Paug);
    xout = zeros(18,36);
    ro = zeros(36,1);
    for j = 1:n
        xout(:,j) = ukf.xaug + sqrt(n+kappa)*S(:,j);
        xout(:,j+n) = ukf.xaug - sqrt(n+kappa)*S(:,j);
        ro(j) = 1/(2*(n+kappa));
        ro(j+n) = 1/(2*(n+kappa));
    end
    ukf.SP = [ukf.xaug xout];
    ukf.ro = [(kappa/(n+kappa)); ro];
end