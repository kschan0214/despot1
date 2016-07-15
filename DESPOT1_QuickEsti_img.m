function [T1map, Momap] = DESPOT1_QuickEsti_img(S1,S2,FA,TR)

FAmap(:,:,1) = ones(size(S1))*FA(1);
FAmap(:,:,2) = ones(size(S2))*FA(2);

img(:,:,1) = S1;
img(:,:,2) = S2;

x = img./tand(FAmap);
y = img./sind(FAmap);

y_diff = y(:,:,2) - y(:,:,1);
x_diff = x(:,:,2) - x(:,:,1);

m = y_diff./x_diff;

T1 = -TR./log(m);
T1(isnan(T1)) = 0;
T1(T1>5) = 5;
T1map = abs(T1);

m_new = exp(-TR./T1);
Momap = mean((y-repmat(m_new,1,1,2).*x)./(1-repmat(m_new,1,1,2)),3);

end