%%  Calculation of DR x f for optimisation of flip angles
%
%   Input:
%           -   Mo: Proton density
%           -   FA: Flip angle
%           -   T1: (same unit as TR)
%           -   TR: (same unit as T1)
%   Output:
%           -   res: Result of (dynamic raneg x signal fraction based on 1 flip angle provided)
%           -   FA1: Same as input FA
%           -   FA2: Paired flip angle provides same signal as FA1
%
%   Author: Kwok-shing Chan @ University of Aberdeen
%   Dtae created: Jan 1, 2016
%   Ref: Rapid combined T1 and T2 mapping using gradient recalled
%   acquisition in the steady state, Deoni et al. MRM 2003;49:515-526
%
function [res, FA1, FA2] = DRxf(Mo, FA, T1, TR)
    [FA1, FA2] = pair_angle(FA, T1, TR);
    DR = normailised_dynamic_range(Mo, FA1, FA2, T1, TR);
    f = singal_fraction(Mo, FA, T1, TR);
    
    res = DR .* f;
end