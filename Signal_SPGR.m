%% Calculation of signal intensity with known flip angle, T1 and TR
%   Input:
%           -   Mo: Proton density, can be arbitary unit
%           -   FA: Flip angle
%           -   T1: (same unit as TR)
%           -   TR: (same unit as T1)
%   Output:
%           -   S: T1-weighted signal
%
%   Author: Kwok-shing Chan @ University of Aberdeen
%   Date created: Jan 1, 2016
%   Ref: Rapid combined T1 and T2 mapping using gradient recalled
%   acquisition in the steady state, Deoni et al. MRM 2003;49:515-526
%
function S = Signal_SPGR(Mo, FA, T1, TR)
E1 = exp(-TR/T1);
S = (Mo.*(1-E1).*sind(FA))./(1-E1.*cosd(FA));
end