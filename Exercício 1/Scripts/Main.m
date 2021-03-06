%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%              - Exercício Computacional de MP208 -             %
%    --- Optimal Filtering with Aerospace Applications ---      %
%                                                               %
%              Autor: João Filipe R. P de A. Silva              %
%                                                               %
%              Main Script: Análise de Estimadores              %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cache Clean-up
clear all
close all
clc

%% General Variables declaration

est.R = 0.01;               %Disturbance variance
est.W = 1;                  %Least Squares Weighing Matrix
est.P_t = [0.01 0;0 0.04];  %Theta Covariance Matrix
est.m_t = [1 2]';           %Theta Expected Value
est.theta = [0 0]';         %Parameter to be estimated
est.N = 100000;              %Number of measurements
est.H = [1:est.N; ones(size(1:est.N))]';

%% Least Square Variables Initialization

ls.A = zeros(2,2);              %Auxiliary LS Sum Matrix
ls.B = zeros(2,1);              %Auxiliary LS Sum Matrix
ls.thetaLS = zeros(2,1);        %Least Squares estimated parameter
ls.RMS_LS = zeros(2,2);         %Least Squares estimator RMS Error

%% Maximum a Posteriori Variables Initialization

map.A = zeros(2,2);             %Auxiliary MAP Sum Matrix
map.B = zeros(2,1);             %Auxiliary MAP Sum Matrix
map.thetaMAP = zeros(2,1);      %Maximum a Posteriori estimated parameter
map.RMS_MAP = zeros(2,2);       %Maximum a Posteriori RMS Error

%% Simulation

for n = 1:100                                          %Monte Carlo Loop
    est.theta = sqrtm(est.P_t)*randn(2,1) + est.m_t;    %Parameter values definition
    
    ls = LS(est,ls,n);                  %Least Squares Estimation
    map = MAP(est,map,n);               %Maximum a Posteriori Estimation
end

%% Sample Means

ls.thetaLS = ls.thetaLS';       %LS Estimative matrix manipulation (to facilitate the mean calculation)
mLS = (sum(ls.thetaLS)/length(ls.thetaLS))'         %LS Estimatives Mean calculation

map.thetaMAP = map.thetaMAP';   %MAP Estimative matrix manipulation (to facilitate the mean calculation)
mMAP = (sum(map.thetaMAP)/length(map.thetaMAP))'    %MAP Estimatives Mean calculation

%% Root Mean Square Error

ls.RMS_LS = sqrtm(ls.RMS_LS/(n-1));      %LS RMS Error calculation
ls.RMS_LS
map.RMS_MAP = sqrtm(map.RMS_MAP/(n-1));  %MAP RMS Error calculation
map.RMS_MAP
sqrtm(map.RMStheo)

%% Histogram Plotting

hLS = histogram2(ls.thetaLS(:,1),ls.thetaLS(:,2));
hLS.BinWidth = [0.1 0.2];
xlabel('\theta_1');
ylabel('\theta_2');
zlabel('Estimate Realizations');

figure

hMAP = histogram2(map.thetaMAP(:,1),map.thetaMAP(:,2));
hMAP.BinWidth = [0.1 0.2];
xlabel('\theta_1');
ylabel('\theta_2');
zlabel('Estimate Realizations');
