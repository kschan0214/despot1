%% Calculation of DR x FS for optimisation of flip angles
%
%   Input:
%           -   Mo: Proton density
%           -   FA1: 1st flip angle
%           -   FA2: 2nd flip angle
%           -   T1: (same unit as TR)
%           -   TR: (same unit as T1)
%   Output:
%           -   res: Result of (dynamic range x fractional signal based on 2 flip angle provided)
%
%   For DR x FS maximisation
%   Author: Kwok-shing Chan @ University of Aberdeen
%   Date created: Jan 1, 2016
%   Ref: Rapid combined T1 and T2 mapping using gradient recalled
%   acquisition in the steady state, Deoni et al. MRM 2003;49:515-526
%
function res = DRxFS(Mo, FA1, FA2, T1, TR)
    DR = normailised_dynamic_range(Mo, FA1, FA2, T1, TR);
    FS = fractional_signal(Mo, FA1, FA2, T1, TR);
    res = DR .* FS;
end







