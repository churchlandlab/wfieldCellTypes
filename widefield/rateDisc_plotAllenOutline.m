function rateDisc_plotAllenOutline(ax,side,allRegions)
%function to add area outlines from allen CCF to current image

load('allenDorsalMapSM.mat')
if ~exist('ax','var') || isempty(ax)
    ax = gca;
end

cIdx = true(1,length(dorsalMaps.edgeOutlineSplitRed));
if exist('side','var')
    if strcmpi(side,'L') || strcmpi(side,'R')
        cIdx = ismember(dorsalMaps.sidesSplitRed,side)';
    end
end

if ~exist('allRegions','var')
    allRegions = false;
end

if ishold(ax)
    checker = true;
else
    hold(ax,'on'); checker = false;
end

%check image size
for x = 1 : length(ax.Children)
    if contains(class(ax.Children(x)),'Image')
        lineScale = min((size(dorsalMaps.allenMask) ./ size(ax.Children(x).CData))); %scale lines to match size of the image
    end
end

for x = find(cIdx)
    if allRegions
        plot(dorsalMaps.edgeOutlineSplit{x}(:,2)./lineScale, dorsalMaps.edgeOutlineSplit{x}(:,1)./lineScale, 'w', 'LineWidth', 0.1);
    else
        plot(dorsalMaps.edgeOutlineSplitRed{x}(:,2)./lineScale, dorsalMaps.edgeOutlineSplitRed{x}(:,1)./lineScale, 'w', 'LineWidth', 0.1);
    end
end

if ~checker
    hold(ax,'off');
end