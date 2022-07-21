%% [t1, m0, mask_fitted] = despot1_mapping(img,fa,tr,mask,b1)
%
% Input
% --------------
% img           : magnitude image, 4D [row,col,slice,flip_angle]
% fa            : flip angle in degree
% tr            : repetition time, either ms or s
% mask          : signal mask, optional
% b1            : B1 map, optional
%
% Output
% --------------
% t1            : T1 map, unit depends on tr
% m0            : proton density map
% mask_fitted   : mask of fitted voxels without encountering error
%
% Description: T1 mapping using DESPOT1 formalism
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 30 Oct 2020
% Date modified: 17 Nov 2020
% Date modified: 21 July 2022
%
function [t1, m0, mask_fitted] = despot1_mapping(img,fa,tr,mask,b1)

dims = size(img);

% if B1 is not provided
if nargin < 5 || isempty(b1)
    b1 = ones(dims(1:3));
end
% if mask is not provided
if nargin < 4 || isempty(mask)
    mask = ones(dims(1:3));
end

% reshape for voxelwise operation
img     = reshape(abs(img),[prod(dims(1:3)) dims(4)]);
mask    = reshape(mask,[prod(dims(1:3)) 1]);
b1      = reshape(b1,[prod(dims(1:3)) 1]);

% masking
ind = find(mask>0);
img = img(ind,:);
b1  = b1(ind);

% DESPOT1
t1 = zeros([prod(dims(1:3)) 1]);
m0 = zeros([prod(dims(1:3)) 1]);
for k = 1:size(img,1)
    if sum(img(k,:)) ~= 0
        [t1(ind(k)), m0(ind(k))] = DESPOT1(img(k,:),fa,tr,'b1',b1(k));
    else
        % avoid matrix inversion with zero matrix
        t1(ind(k)) = 0; m0(ind(k)) = 0;
    end
end
t1      = reshape(t1,dims(1:3));
m0      = reshape(m0,dims(1:3));
mask    = reshape(mask,dims(1:3));

% remove physically infeasible values
mask_fitted            = ones(size(t1));
mask_fitted(t1<=0)     = 0;                 % T1=0 is also physically not feasible
mask_fitted(isnan(t1)) = 0;
mask_fitted(isinf(t1)) = 0;
mask_fitted(m0<0)      = 0;
mask_fitted(isnan(m0)) = 0;
mask_fitted(isinf(m0)) = 0;
mask_fitted            = mask_fitted .* mask;

t1 = t1 .* mask_fitted;
m0 = m0 .* mask_fitted;
 
% t1(t1<0) = 0;
% t1(isnan(t1)) = 0;
% t1(isinf(t1)) = 0;
% 
% m0(m0<0) = 0;
% m0(isnan(m0)) = 0;
% m0(isinf(m0)) = 0;

end