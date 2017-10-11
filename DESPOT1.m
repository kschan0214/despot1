%% [T1, Mo] = DESPOT1(S, FA, TR,varargin)
%
% Usage:
%   [T1, Mo] = DESPOT1(S, FA, TR,'range',rangeT1,'option',options);
% 
% Input
% ----------
% S         : Set of T1-weighted signals from different flip angles, pixel-based
% FA        : Set of flip angles used corresponding to S
% TR        : TR of sequence (s)
% Flags     :
%   'range'     -   boundary of T1
%   'option'    -   fitting option
%
% Output
% ----------
% T1        : Calcuted T1 (same as TR)
% Mo        : Calculated proton density
%
% DESPOT1 to estimate T1 and proton density Mo
% Linear fitting of y = ax+b
% Rapid combined T1 and T2 mapping using gradient recalled acquisition in the steady state.
% Deoni et al. MRM 2003;49:515-526
%
% Kwok-shing Chan @ dccn
% k.chan@donders.ru.nl
% Date created: Jan 11, 2016
% Date last edited: July 14, 2017
%
function [T1, Mo] = DESPOT1(S, FA, TR,varargin)
    if length(S) < 2
        error('Only one point to fit.');
    end
    % parse arhument input
    [T1_lb,T1_ub,option] = parse_varargin_DESPOT1(varargin);
    
    %% Obtain initial guesses
    [T10, Mo0] = DESPOT1_QuickEsti(abs(S),FA,TR,T1_ub,T1_lb);
    c0 = [T10, Mo0];
    lb = [T1_lb, min(S)];
    ub = [T1_ub, 2*Mo0];
    
    res = lsqnonlin(@(x)fitError_DESPOT1(x,S,FA,TR), c0,lb,ub,option);
    T1 = abs(res(1));
    Mo = abs(res(2));
end

%% fast estimation by fitting to straightline
function [T1map, Momap] = DESPOT1_QuickEsti(S,FA,TR,T1_ub,T1_lb)

x = S./tand(FA);
y = S./sind(FA);

y_diff = y(2) - y(1);
x_diff = x(2) - x(1);

m = y_diff./x_diff;

T1 = -TR./log(m);
T1(isnan(T1)) = T1_lb;
T1(isinf(T1)) = T1_ub;
T1(T1<T1_lb) = T1_lb;
T1(T1>T1_ub) = T1_ub;
T1map = abs(T1);

m_new = exp(-TR./T1);
Momap = mean((y-repmat(m_new,2,1).*x)./(1-repmat(m_new,2,1)));
Momap(isnan(Momap)) = 0;
Momap(isinf(Momap)) = 0;
end

%% compute fitting residual
function [fiter] = fitError_DESPOT1(x,S_meas,FA, TR)
T1 = x(1);
Mo = x(2);

E1 = exp(-TR/T1);
S_fit = E1.*(S_meas./tand(FA))+Mo*(1-E1);
fiter = abs(S_fit - S_meas./sind(FA));
end

%% parse argument input
function [T1_lb,T1_ub,option] = parse_varargin_DESPOT1(arg)
T1_lb = 0;
T1_ub = 5;
fitOption = false;
if ~isempty(arg)
    for kvar = 1:length(arg)
        if strcmpi(arg{kvar},'range')
            range = arg{kvar+1};
            T1_lb = range(1);
            T1_ub = range(2);
        end
        if strcmpi(arg{kvar},'option')
            option = arg{kvar+1};
            fitOption = true;
        end
    end
    if fitOption == false
        option = optimoptions(@lsqnonlin,'Display','off','Jacobian','off','DerivativeCheck','off','MaxIter',1000);
    end
end
end
