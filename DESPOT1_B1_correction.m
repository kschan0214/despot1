%% DESPOT1 to estimate T1 and proton density Mo (unfinished)
%   
%   Linear fitting of y = ax+b
%   Rapid combined T1 and T2 mapping using gradient recalled acquisition in the steady state.
%   Deoni et al. MRM 2003;49:515-526
%
%   Creator: Kwok-shing Chan @ University of Aberdeen
%   Date created: Jan 11, 2016
%
%%
function outParams = DESPOT1_B1_correction(imDataParams)
im = double(imDataParams.images);
[sx,sy,nFA] = size(im);

FA = imDataParams.FA;
TR = imDataParams.TR;
flag = imDataParams.flag;

% Get the signal norm at each voxel (used for thresholding)
smagn = sqrt(sum(sum(abs(im).^2,5),4));
smagn = smagn/max(smagn(:));
THRESHOLD = 0.05;

% Initialize maps and iterate over voxels
T1Map = zeros(sx,sy);
B1Map = zeros(sx,sy);
resnorm = zeros(sx,sy);

for kx=1:sx
    for ky=1:sy
        if mod(kx,10) ==0
            display(kx);
        end
        if smagn(kx,ky) > THRESHOLD
            [c, resnorm(kx, ky)] = DESPOT1_B1corr(squeeze(im(kx,ky,:)),FA, TR, flag(kx, ky));
            T1Map(kx,ky) = c(2);
            B1Map(kx,ky) = c(1);
        end
    end
end
outParams.T1Map = T1Map;
outParams.B1Map = B1Map;

function [c, resnorm] = DESPOT1_B1corr(S, FA, TR, flag)
    if length(S) < 2
        display('Only one point to fit.');
    end
    if nargin <4
        flag = 0;
    end
    
    opts = optimset('lsqcurvefit');
    options = optimset(opts,'Display','off','MaxIter',10);
    if flag == 0
    %% Test for fat T1
    c0 = [1, 371];
    lb = [0, 250];
    ub = [2, 420];
    [c1, resnorm1] = lsqnonlin(@(x)DESPOT1_B1correction(x, S, FA, TR),c0,lb,ub,options);
    %% Test for water T1
    c0 = [1, 1420];
    lb = [0, 1000];
    ub = [2, 1500];
    
    [c2, resnorm2] = lsqnonlin(@(x)DESPOT1_B1correction(x, S, FA, TR),c0,lb,ub,options);
    
    %% Comparing the sos residue
    if resnorm1 <= resnorm2
        c = c1;
        resnorm = resnorm1;
    else 
        c = c2;
        resnorm = resnorm2;
    end
    elseif flag == 2
    c0 = [1, 365];
    lb = [0, 360];
    ub = [2, 370];
    [c, resnorm] = lsqnonlin(@(x)DESPOT1_B1correction(x, S, FA, TR),c0,lb,ub,options);
    
    elseif flag == 1
    c0 = [1, 1420];
    lb = [0, 1350];
    ub = [2, 1500];
    [c, resnorm] = lsqnonlin(@(x)DESPOT1_B1correction(x, S, FA, TR),c0,lb,ub,options);
    end

function F = DESPOT1_B1correction(x, S, FA, TR)
B1 = x(1);
T1 = x(2);
E1 = exp(-TR/T1);

I = (sind(B1.*FA).*(1-E1))./(1-E1.*cosd(B1.*FA));
I_mat = repmat(I,1,length(FA));
I_inv = diag(1./I);
res_est = I_mat*I_inv;

S_mat = repmat(S,1,length(S));
S_inv = diag(1./S);
res_meas = S_mat*S_inv;

F = reshape(abs(res_est-res_meas),[length(res_meas(:)), 1]);

