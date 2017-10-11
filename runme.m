%% This script demostrates how to use DESPOT1 for T1 mapping

addpath func/
addpath utils/

%% Simulate signal
M0 = 100;       % a.u.
FA = [5,30,60]; % degree
T1 = 800e-3;    % s
TR = 50e-3;     % s
S = Signal_GRE_T1wMono(M0,FA,T1,TR);

%% using DESPOT1
% regression method
[t1,m0] = DESPOT1(S,FA,TR);
% lsqnonneg method
[t1,m0] = DESPOT1(S,FA,TR,'lsqnonneg');
% lsqnonlin method
option = optimoptions(@lsqnonlin,'Display','off','Jacobian','off','DerivativeCheck','off','MaxIter',100);
[t1,m0] = DESPOT1(S,FA,TR,'lsqnonlin','option',option);
