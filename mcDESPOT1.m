%% function output = function_name(input)
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
function [res,resnorm]=mcDESPOT1(S,FA,TR,varargin)
% parse argument input
[initGuess,option] = parse_varargin_mcDESPOT1(varargin);

%% Core
% DESPOT1 formulation
S = double(S);
FA = double(FA);
TR = double(TR);
    
[t10,rho0] = DESPOT1(S,FA,TR);

% solve DESPOT1 using different approach
if isempty(initGuess)
    c0=[100e-3,rho0*0.2,t10,rho0*0.8];
else
    c0=initGuess(:).';
end
lb = [50e-3,0,250e-3,0];
ub = [250e-3,2*rho0,800e-3,rho0*2];

        [estimate,resnorm] = lsqnonlin(@(x)fitError_mcDESPOT1(x,S,FA,TR),c0,lb,ub,option);
        res(1).t1 = estimate(1);
        res(1).rho = estimate(2);
        res(2).t1 = estimate(3);
        res(2).rho = estimate(4);

end

%% lsqnonlin: compute fitting residual 
function [fiter] = fitError_mcDESPOT1(x,S_meas,FA,TR)
% grab estimates
t1_a = x(1);
rho_a = x(2);
t1_b = x(3);
rho_b = x(4);

% simulate signal
S_fit = Signal_GRE_T1wMono(rho_a,FA,t1_a,TR) + Signal_GRE_T1wMono(rho_b,FA,t1_b,TR) ;
% compute fiter, using magnitude fitting
fiter = computeFiter(S_meas,S_fit,length(S_fit));

end


%% parse argument input
function [initGuess,option] = parse_varargin_mcDESPOT1(arg)
% predefine parameters
option = [];
initGuess = [];

% look for flags and 'Name/Value' pairs
for kvar = 1:length(arg)
    if strcmpi(arg{kvar},'init')
        initGuess = arg{kvar+1};
    end
    if strcmpi(arg{kvar},'option')
        option = arg{kvar+1};
    end
end
% lsqnonlin fitting option
if isempty(option)
    option = optimoptions(@lsqnonlin,'Display','off','Jacobian','off','DerivativeCheck','off','MaxIter',500);
end
end