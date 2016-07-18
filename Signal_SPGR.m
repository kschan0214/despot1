%% Calculation of signal intensity with known flip angle, T1 and TR
%   S = Signal_SPGR(Mo, FA, T1, TR, necho, TE0, dTE, T2s)
%   Theoretical T1-weighted GRE signal
%   Compaticale with multi-echo signal acquisition
%
%   Input:
%           -   Mo: Proton density, can be arbitary unit
%           -   FA: Flip angle
%           -   T1: (same unit as TR)
%           -   TR: (same unit as T1)
%           -   necho: number of echo(es) (same unit as T1) (optional)
%           -   TE0: initial TE (same unit as T1) (optional)
%           -   dTE: echo spacing (same unit as T1) (optional)
%           -   T2s: T2* of sample (same unit as T1) (optional)
%   Output:
%           -   S: T1-weighted signal
%
%   Author: Kwok-shing Chan @ University of Aberdeen
%   Date created: Jan 1, 2016
%   Ref: Rapid combined T1 and T2 mapping using gradient recalled
%   acquisition in the steady state, Deoni et al. MRM 2003;49:515-526
%   Date last edited: July 18, 2016
%
function S = Signal_SPGR(Mo, FA, T1, TR, necho, TE0, dTE, T2s)
if nargin<5
    TE0 = 0;
    dTE = 0;
    necho = 0;
    T2s = eps;
end
if necho ~= 0
    te = [TE0];
    for n=1:necho-1
        te = [te;n*dTE+TE0];
    end
else
    te = 0;
end
E1 = exp(-TR/T1);
S = (Mo.*(1-E1).*sind(FA))./(1-E1.*cosd(FA)).* exp(-te/T2s);
end