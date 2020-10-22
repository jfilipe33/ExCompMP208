%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%              - Exercício Computacional de MP208 -             %
%    --- Optimal Filtering with Aerospace Applications ---      %
%                                                               %
%              Autor: João Filipe R. P de A. Silva              %
%                                                               %
%        Function Script: Estimador Maximum a Posteriori        %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function map = MAP(sys,map,i)
       
    map.A = map.A + sys.H'*inv(sys.R)*sys.H;
    Pn = inv(map.A + inv(sys.P_t));
    map.B = map.B + sys.H'*inv(sys.R)*sys.y;
    map.thetaMAP(:,i) = Pn*inv(sys.P_t)*sys.m_t + Pn*(map.B);
    map.RMS_MAP = map.RMS_MAP + ((map.thetaMAP(:,i)-sys.theta)*(map.thetaMAP(:,i)-sys.theta)');

    P1 = Pn*map.A - eye(2);
    map.RMStheo = Pn*inv(sys.P_t)*sys.m_t*sys.m_t'*inv(sys.P_t)*Pn + P1*(sys.P_t + sys.m_t*sys.m_t')*P1 + Pn*inv(sys.P_t)*sys.m_t*sys.m_t'*P1 + P1*sys.m_t*sys.m_t'*inv(sys.P_t)*Pn + Pn*map.A*Pn; 
end