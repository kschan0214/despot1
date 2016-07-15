%% Calculation of signal fraction when S1 may not equal to S2
%
%   Input:
%           -   Mo: Proton density
%           -   FA1: 1st flip angle
%           -   FA2: 2nd flip angle
%           -   T1: (same unit as TR)
%           -   TR: (same unit as T1)
%   Output:
%           -   FS: Signal fraction of using FA1 and FA2 compared to Ernst
%                   angle
%
%   Author: Kwok-shing Chan @ University of Aberdeen
%   Date created: Jan 1, 2016
%   Ref: Rapid combined T1 and T2 mapping using gradient recalled
%   acquisition in the steady state, Deoni et al. MRM 2003;49:515-526
%
function FS = fractional_signal(Mo, FA1, FA2, T1, TR)
    S1 = Signal_SPGR(Mo, FA1, T1, TR);
    S2 = Signal_SPGR(Mo, FA2, T1, TR);
    EA = ernst_angle(T1,TR);
    SE = Signal_SPGR(Mo, EA, T1, TR);
    FS = (S1 + S2) ./ (2.*SE);
end