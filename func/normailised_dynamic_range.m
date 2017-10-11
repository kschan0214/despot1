%% Calculation of normalised dynamic range of regression line 
%
%   Input:
%           -   Mo: Proton density
%           -   FA1: 1st flip angle
%           -   FA2: 2nd flip angle
%           -   T1: (same unit as TR)
%           -   TR: (same unit as T1)
%   Output:
%           -   DR: Dynamic range
%
%   Author: Kwok-shing Chan @ University of Aberdeen
%   Date created: Jan 1, 2016
%   Ref: Rapid combined T1 and T2 mapping using gradient recalled
%   acquisition in the steady state, Deoni et al. MRM 2003;49:515-526
%
function DR = normailised_dynamic_range(Mo, FA1, FA2, T1, TR)
    S1 = Signal_SPGR(Mo, FA1, T1, TR);
    S2 = Signal_SPGR(Mo, FA2, T1, TR);
    
    DR = S2./(Mo.*sind(FA2)) - S1./(Mo.*sind(FA1));
end