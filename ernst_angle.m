%% Calculation of Ernst angle with knwon T1 and TR
%% function EA = ernst_angle(T1,TR)
%   Input:
%           -   T1: (same unit as TR)
%           -   TR: (same unit as T1)
%   Output:
%           -   EA: Ernst angle based on T1 and TR
%
%   Author: Kwok-shing Chan @ University of Aberdeen
%   Date created: Jan 1, 2016
%   Ref: Rapid combined T1 and T2 mapping using gradient recalled
%   acquisition in the steady state, Deoni et al. MRM 2003;49:515-526
%
function EA = ernst_angle(T1,TR)
    EA = acosd(exp(-TR/T1));
end