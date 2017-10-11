%%  Calculation of pair flip angle of S1=S2 for given flip angle
%   Input:
%           -   FA: target flip angle
%           -   T1: (same unit as TR)
%           -   TR: (same unit as T1)
%   Output:
%           -   [FA1, FA2]: paired angles based on FA
%
%   Author: Kwok-shing Chan @ University of Aberdeen
%   Date created: Jan 1, 2016
%   Ref: Rapid combined T1 and T2 mapping using gradient recalled
%   acquisition in the steady state, Deoni et al. MRM 2003;49:515-526
%
function [FA1, FA2] = pair_angle(FA, T1, TR)
    E1 = exp(-TR/T1);
    FA1 = zeros(1,length(FA));
    FA2 = zeros(1,length(FA));
    for n=1:length(FA)        
        syms x
        eqn = sin(x*pi/180)/(1-E1*cos(x*pi/180)) == sin(FA(n)*pi/180)/(1-E1*cos(FA(n)*pi/180));
        solx = solve(eqn,x);
        % Put FA into FA1, estimated flip angle into FA2
        if ((FA(n)-abs(double(solx(1))))^2 < (FA(n) - abs(double(solx(2))))^2)
            FA1(n) = abs(double(solx(1)));
            FA2(n) = abs(double(solx(2)));
        else
            FA1(n) = abs(double(solx(2)));
            FA2(n) = abs(double(solx(1)));
        end
    end
end