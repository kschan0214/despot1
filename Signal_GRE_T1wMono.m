%% S = Signal_GRE_T1wMono(Mo, FA, T1, TR,varargin)
%
% Input:
% ------
%	Mo: Proton density, can be arbitary unit
%	FA: Flip angle
% 	T1: (same unit as TR)
%	TR: (same unit as T1)
% Output:
% -------
%	S: T1-weighted signal
%
% Description: Calculation of signal intensity with known flip angle, T1 and TR
%
%   Author: Kwok-shing Chan @ University of Aberdeen
%   Date created: Jan 1, 2016
%   Ref: Rapid combined T1 and T2 mapping using gradient recalled
%   acquisition in the steady state, Deoni et al. MRM 2003;49:515-526
%   Date last edited: July 18, 2016
%
function S = Signal_GRE_T1wMono(Mo,FA,T1,TR,varargin)

if length(T1)>1 && length(FA)>1
    disp('More than one T1. Output will be in 2D.');
end

%% Core algorithm
S = zeros(length(T1),length(FA));
for k=1:length(T1)
    E1 = exp(-TR/T1(k));
    S(k,:) = (Mo.*(1-E1).*sind(FA))./(1-E1.*cosd(FA));
end

%% for simulation display
if ~isempty(varargin)
    if strcmpi(varargin{1,1},'image')
        figure;imagesc(S);xlabel('Flip angle (degree)');ylabel('T1');colormap jet;
        set(gca,'Ydir','normal');
        yTickLabelName = get(gca,'yticklabel');
        for kk = 1:length(yTickLabelName)
            yTickLabelName{kk,1} = num2str(T1(str2double(yTickLabelName{kk,1})));
        end
        set(gca,'yTickLabel',yTickLabelName);
        xTickLabelName = get(gca,'xticklabel');
        for kk = 1:length(xTickLabelName)
            xTickLabelName{kk,1} = num2str(FA(str2double(xTickLabelName{kk,1})));
        end
        set(gca,'xTickLabel',xTickLabelName);
        title(sprintf('SPGR T1 weighting, TR=%f',TR));
    end
    if strcmpi(varargin{1,1},'surf')
        [rx,ry] = meshgrid(FA,T1);
        figure;surf(rx,ry,S);xlabel('Flip angle (degree)');ylabel('T1');colormap jet;
        title(sprintf('SPGR T1 weighting, TR=%f',TR));
    end
    if strcmpi(varargin{1,1},'plot')
        figure;
        for k = 1:length(T1)
            plot(FA,squeeze(S(k,:)));hold on;
        end
        xlabel('Flip angle (degree)');ylabel('Signal intensity (a.u.)');colormap jet;
        title(sprintf('SPGR T1 weighting, TR=%f',TR));
    end
end
end