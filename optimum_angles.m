%% Calculation of optimum flip angles based on S1=S2
%% FA = optimum_angles(T1, TR)
%
%   Input:
%           -   T1: (same unit as TR)
%           -   TR: (same unit as T1)
%   Output:
%           -   FA: optimum flip angles (same unit as FA)
%
%   Author: Kwok-shing Chan @ University of Aberdeen
%   Date created: Jan 1, 2016
%   Ref: Rapid combined T1 and T2 mapping using gradient recalled
%   acquisition in the steady state, Deoni et al. MRM 2003;49:515-526
function FA = optimum_angles(T1, TR)
    f= 0.71;
    E1 = exp(-TR/T1);
    sol1 = (f.^2*E1 + (1-E1.^2).*sqrt(1-f.^2))./(1-E1.^2*(1-f.^2));
    sol2 = (f.^2*E1 - (1-E1.^2).*sqrt(1-f.^2))./(1-E1.^2*(1-f.^2));
    
    FA(1) = acosd(sol1);
    FA(2) = acosd(sol2);
end