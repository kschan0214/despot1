%% [t1,m0]=DESPOT1_regression(S,FA,TR)
%
% Usage:
%
% Input
% --------------
%
% Output
% --------------
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 
% Date last modified:
%
%
function [t1,m0]=DESPOT1_regression(S,FA,TR)
y = S(:)./sind(FA(:));
x = S(:)./tand(FA(:));
B = ones(length(S),2);
B(:,1) = x;

b1 = pinv(B)*y;

t1 = -TR/log(b1(1));

m0 = b1(2)/(1-exp(-TR/t1));

end