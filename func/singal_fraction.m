%% Calculation of signal fraction when S1=S2
%   Input:
%           -   Mo: Proton density, can be arbitary unit
%           -   FA: Flip angle
%           -   T1: (same unit as TR)
%           -   TR: (same unit as T1)
%   Output:
%           -   f: Signal fraction compared to maximum signal (Ernst
%                  angle's signal)
%
%   Author: Kwok-shing Chan @ University of Aberdeen
%   Date created: Jan 1, 2016
%   Ref: Rapid combined T1 and T2 mapping using gradient recalled
%   acquisition in the steady state, Deoni et al. MRM 2003;49:515-526
%
function f = singal_fraction(Mo, FA, T1, TR)
    EA = ernst_angle(T1, TR);
    SE = Signal_SPGR(Mo, EA, T1, TR);
    S = Signal_SPGR(Mo, FA, T1, TR);
    
    f = S./SE; 
end