%% DESPOT1 to estimate T1 and proton density Mo
%
%   Input:
%           -   S: T1-weighted signal
%           -   FA: Set of flip angles used corresponding to S
%           -   TR: TR of protocol (ms)
%           -   flag: Prior knowledge of fat/water distribution 
%                       - 0: No prior knowledge
%                       - 1: Water dominated pixel
%                       - 2: Fat dominated pixel
%   Output:
%           -   T1: Calcuted T1 (ms)
%           -   Mo: Calculated proton density
%
%   Linear fitting of y = ax+b
%   Rapid combined T1 and T2 mapping using gradient recalled acquisition in the steady state.
%   Deoni et al. MRM 2003;49:515-526
%
%   Author: Kwok-shing Chan @ University of Aberdeen
%   Date created: Jan 11, 2016
%
function [T1, Mo] = DESPOT1(S, FA, TR, flag)
    if length(S) < 2
        display('Only one point to fit.');
    end
    if nargin <4
        flag = 0;
    end
    fun = @(c) c(1).*(S./tand(FA)) + c(2) - (S./sind(FA));
    
    opts = optimset('lsqcurvefit');
    options = optimset(opts,'Display','off','MaxIter',20);
    if flag == 0
    %% Test for fat T1
    E1f_0 = exp(-TR/371);   E1_lb = exp(-TR/80);  E1_ub = exp(-TR/2000);
    Mof_0 = mean(S.*((1-E1f_0.*cosd(FA))./((1-E1f_0).*sind(FA))));

    c0 = [E1f_0, (1-E1f_0)*Mof_0];
    lb = [E1_lb, (1-E1_ub)*Mof_0];
    ub = [E1_ub, (1-E1_lb)*Mof_0];
    [c1, resnorm1] = lsqnonlin(fun,c0,lb,ub,options);
    %% Test for water T1
    E1w_0 = exp(-TR/1275);
    Mow_0 = mean(S.*((1-E1w_0.*cosd(FA))./((1-E1w_0).*sind(FA))));

    c0 = [E1w_0, (1-E1w_0)*Mow_0];
    lb = [E1_lb, (1-E1_ub)*Mow_0];
    ub = [E1_ub, (1-E1_lb)*Mow_0];
    
    [c2, resnorm2] = lsqnonlin(fun,c0,lb,ub,options);
    
    %% Comparing the sos residue
    if resnorm1 <= resnorm2
        c = c1;
    else 
        c = c2;
    end
    elseif flag == 2
        E1f_0 = exp(-TR/371);   E1_lb = exp(-TR/80);  E1_ub = exp(-TR/500);
        Mof_0 = mean(S.*((1-E1f_0.*cosd(FA))./((1-E1f_0).*sind(FA))));

        c0 = [E1f_0, (1-E1f_0)*Mof_0];
        lb = [E1_lb, (1-E1_ub)*Mof_0];
        ub = [E1_ub, (1-E1_lb)*Mof_0];
        [c] = lsqnonlin(fun,c0,lb,ub,options);
    elseif flag == 1
        E1w_0 = exp(-TR/1275); E1_lb = exp(-TR/500);  E1_ub = exp(-TR/2000);
        Mow_0 = mean(S.*((1-E1w_0.*cosd(FA))./((1-E1w_0).*sind(FA))));

        c0 = [E1w_0, (1-E1w_0)*Mow_0];
        lb = [E1_lb, (1-E1_ub)*Mow_0];
        ub = [E1_ub, (1-E1_lb)*Mow_0];
    
        [c] = lsqnonlin(fun,c0,lb,ub,options);
    end
    m = c(1);
    b = c(2);
    T1 = -TR/log(m);
    Mo = b/(1-m);
end