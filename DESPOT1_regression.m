%% [t1,m0]=DESPOT1_regression(S,FA,TR)
%
% Input
% -----
% S     : 1D signal
% FA    : 1D flip angle
% TR    : repetition time
%
% Output
% ------
% t1    : T1
% m0    : proton density weighted signal
%
% Description: Solving DESPOT1 equation by linear regression
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 11 October 2017
% Date last modified:
%
%
function [t1,m0]=DESPOT1_regression(S,FA,TR)
S = double(S);
FA = double(FA);
TR = double(TR);

y = S(:)./sind(FA(:));
xCol = S(:)./tand(FA(:));
x = ones(length(S),2);
x(:,1) = xCol;

b = x\y;
% b = lsqnonneg(x,y);

t1 = -TR/log(b(1));

m0 = b(2)/(1-exp(-TR/t1));

end