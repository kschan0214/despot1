%% Calculation of optimum flip angles with known T1 and TR
%
%   Input:
%           -   T1: (same unit as TR)
%           -   TR: (same unit as T1)
%           -   range: range of flip angle allowed (optional)
%   Output:
%           -   FA_opt: optimum flip angles (same unit as FA)
%
%   Author: Kwok-shing Chan @ University of Aberdeen
%   Date created: Jan 1, 2016
%   Ref: Rapid combined T1 and T2 mapping using gradient recalled
%   acquisition in the steady state, Deoni et al. MRM 2003;49:515-526
%
function FA_opt = optimum_angles_brust(T1, TR, range)
% (Optional)Calculate the optimum flip angle within the possible range of
% flip angle
if nargin < 3
    lb = 0.1;
    resolution = 0.05;
    ub = 180;
else
    lb = range(1);
    resolution = 0.05;
    ub = range(2);
end
Mo = 1e6;
[FA1, FA2] = ndgrid(lb:resolution:ub, lb:resolution:ub);

res = DRxFS(Mo, FA1, FA2, T1, TR);
[m, inx] = max(res(:));
[y, x] = ind2sub(size(res),inx);
FA_opt(1) = FA1(y,x);
FA_opt(2) = FA2(y,x);
FA_opt = sort(FA_opt);
end