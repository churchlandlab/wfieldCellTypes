% code to analyze umap clustering of single-region locaNMF components
load('allenDorsalMapSM.mat')
allenMask = dorsalMaps.allenMask;
dataOverview = rateDiscRecordings;
groups = {'mSM' 'Fez' 'Plex' 'CSP'};
nrNeighbors = 20; %nr of neigbours for classifier
pcCnt = 1000; %maximum nr of components from each animal per run
nrRuns = 100; %nr of repetitions for classifier performance
noCStr = true;
groupColor = {[1 0 0] [0 0 1] [0 1 0] [0 1 1]}; %colors for groups
Markers = {'+','o','*','x','v','d','^','s','>','<', '+','o','*'};
rng(0);

cPath = 'X:\RawData\'; %path to raw data
umapPath = 'X:\RawData\Celltypes_dataset\umap\'; %path to umap results
umapFilename = '20220529T1912_localA'; %base name of umap files

%% show locanMF results
load([cPath 'localRegionMap.mat'])
regionMap = localRegionMap;
[~, regionLabels] = rateDisc_makeRegionMap; %get regiolabels
lHS = ~allenMask;
lHS(:, 1: round(size(allenMask,2)/2)+1) = false;
regionMapPlot = regionMap;
regionMapPlot = arrayShrink(regionMapPlot,allenMask,'merge');
regionMapPlot = arrayShrink(regionMapPlot,allenMask,'split');
regionMapPlot(regionMapPlot == 0) = NaN;
regionMapPlot(lHS) = regionMapPlot(lHS)+max(regionMapPlot(:));
for iRegions = (1:10) + max(regionMap(:))
    cFrame = regionMapPlot == iRegions;
    cFrame = imclose(cFrame,strel('disk',6));
    regionMapPlot(cFrame) = iRegions;
end

figure('renderer','painters');
% subplot(1,4,1);
cImg = imageScale(regionMapPlot);
rateDisc_plotAllenOutline(gca,'L');
colormap(cImg.Parent,hsv(256));
caxis([0 max(regionMapPlot(:))]); hold on;
for iRegions = (1:max(regionMap(:))) + max(regionMap(:))
    cFrame = regionMapPlot == iRegions;
    a = bwboundaries(cFrame); %outline of selected area
    areaOutline{iRegions} = a{1};
    for x = 1 : length(areaOutline)
        plot(smooth(areaOutline{iRegions}(:,2),10),smooth(areaOutline{iRegions}(:,1),10),'k', 'linewidth', 2)
    end
end

%% show locaNMF examples
load([cPath 'Fez71\SpatialDisc\16-Jul-2020\newAC_20_50.mat'], 'A')
figure('renderer','painters');
subplot(3,1,1);
fezExp = A(:,:,76);
fezExp(isnan(fezExp)& ~allenMask) = 0;
imageScale(arrayCrop(fezExp, allenMask)); colormap inferno; axis image
title('Fez - Fez71 - 16-Jul-2020 - 4'); rateDisc_plotAllenOutline(gca);
caxis([0 0.75])
 
subplot(3,1,2);
fezExp = A(:,:,104);
fezExp(isnan(fezExp)& ~allenMask) = 0;
imageScale(arrayCrop(fezExp, allenMask)); colormap inferno; axis image
title('Fez - Fez71 - 16-Jul-2020 - 4'); rateDisc_plotAllenOutline(gca);
caxis([0 0.75])

subplot(3,1,3);
fezExp = A(:,:,115);
fezExp(isnan(fezExp)& ~allenMask) = 0;
imageScale(arrayCrop(fezExp, allenMask)); colormap inferno; axis image
title('Fez - Fez71 - 16-Jul-2020 - 4'); rateDisc_plotAllenOutline(gca);
caxis([0 0.75])

%% get umap results from all regions
if noCStr
    rmvGroup = 4; %this is the group that should not be considered (4 for CSP neurons)
else
    rmvGroup = 0;
end

% new version with temporal components
umapFile = dir([umapPath '*' umapFilename '_umap_allAreas.mat']); % load this for local - locaNMF components used in first submission
tempUmapFile = dir([umapPath '*' umapFilename '_umap_allAreas_temporal.mat']); % load this for local - locaNMF components
compFile = dir([umapPath '*' umapFilename '.mat']); % load this for local - locaNMF components

load([umapPath tempUmapFile.name], 'umapOut')
umapOut(:,1) = umapOut(:,1) - min(umapOut(:,1));
umapOut(:,2) = umapOut(:,2) - min(umapOut(:,2));
umapOut = round(umapOut * 100);
tempUmapOut = umapOut;

load([umapPath umapFile.name], 'umapOut')
load([umapPath compFile.name], 'groupIdx', 'mouseIdx', 'areaIdx', 'X', 'Y')
Y = Y';

rmvIdx = groupIdx == rmvGroup;
umapOut(:,1) = umapOut(:,1) - min(umapOut(:,1));
umapOut(:,2) = umapOut(:,2) - min(umapOut(:,2));
umapOut = round(umapOut * 100);

X = X(:,~rmvIdx);
Y = Y(:,~rmvIdx);
umapOut = umapOut(~rmvIdx,:);
tempUmapOut = tempUmapOut(~rmvIdx,:);
mouseIdx = mouseIdx(~rmvIdx);
groupIdx = groupIdx(~rmvIdx);
areaIdx = areaIdx(~rmvIdx);

%% show umap projection for all regions
h1 = figure('renderer','painters'); hold on;
cMice = unique(mouseIdx);
compCnt = size(umapOut,1); %nr of components
stepSize = ceil(compCnt / 100); %determines how many PCs are processed at a time
xx = 1 : stepSize : compCnt;

% need to do some shuffling here to not get an artifical clustering pattern
% by overlaying groups on each other
for shuffleRuns = xx(randperm(length(xx)))

    stepIdx = [shuffleRuns, shuffleRuns + stepSize];
    if stepIdx(end) > compCnt; stepIdx(end) = compCnt; end
    
    shuffIdx = randperm(length(cMice)); %shuffle order of mice to avoid certain groups to dominate

    for x = shuffIdx
        
        useIdx = mouseIdx == cMice(x);
        cGroup = unique(groupIdx(useIdx));
        
        useIdx([1:stepIdx(1) stepIdx(end):end]) = false;
        
        if sum(useIdx) > 0
            
            scatter(umapOut(useIdx,1), umapOut(useIdx,2), [], groupIdx(useIdx), ...
                    'Marker', Markers{x}, 'SizeData', 20, 'MarkerEdgeColor', groupColor{cGroup}); axis square;
            
        end
    end
end


%% new cell type classifier - based on separate umap projections
baseName = strrep(compFile.name, '.mat', '_mouseTest'); %base name for umap files
load([umapPath [baseName '1.0.mat']], 'mouseIdx', 'groupIdx')

allGroups = unique(groupIdx(groupIdx ~= rmvGroup));
nrGroups = length(allGroups);
cMice = unique(mouseIdx(groupIdx ~= rmvGroup));
estCompClass = cell(1,length(cMice));
compRecIdx = cell(1,length(cMice));
decoderPerf = cell(1, nrGroups);
for iMice = 1 : length(cMice)
        
    cMouse = cMice(iMice);
    cFile = [baseName num2str(cMouse) '.0.mat'];
    load([umapPath cFile], 'umapRed', 'umapTest', 'recIdx' , 'mouseIdx', 'areaIdx', 'groupIdx')

    compRecIdx{iMice} = recIdx(mouseIdx == cMouse); %keep IDs for each moouse to make sure you can find them again
    cGroup = unique(groupIdx((mouseIdx == cMouse)));
    testPCs = size(umapTest, 1);

    groupIdx = groupIdx(mouseIdx ~= cMouse);
    areaIdx = single(areaIdx(mouseIdx ~= cMouse));
    useIdx = groupIdx ~= rmvGroup; %dont use rmvGroup
    allGroups = unique(groupIdx(useIdx));
    
    %get index with even amount of regions from each cell types
    evenRegionIdx = false(size(areaIdx));
    minCnt = inf;
    for iGroups = allGroups'
        [~, cCnt] = rateDisc_animalIndex(areaIdx(groupIdx == iGroups & useIdx), inf); %get even nr of components from all regions
        minCnt = min([minCnt, cCnt]);
    end
    for iGroups = allGroups'
        cAreaIdx = areaIdx;
        cAreaIdx(groupIdx ~= iGroups | ~useIdx) = NaN;
        cAreaIdx = rateDisc_animalIndex(cAreaIdx, minCnt); %get even nr of components from all regions
        evenRegionIdx(cAreaIdx) = true;
    end
    
    % compute performance
    othGroupIdx = groupIdx;
    othGroupIdx(~evenRegionIdx | ~useIdx) = NaN;
    estCompClass{iMice} = NaN(nrRuns, testPCs, 'single');
    for iRuns = 1 : nrRuns
        
        [othIdx, minCnt] = rateDisc_animalIndex(othGroupIdx, pcCnt); %get even nr of components from all groups
        othUmap = umapRed(othIdx, :);
        othGroups = othGroupIdx(othIdx);

        cID = NaN(nrNeighbors, testPCs);
        for iComps = 1 : testPCs
            a = umapTest(iComps, :);
            cDist = sqrt((a(1)-othUmap(:,1)).^2+(a(2)-othUmap(:,2)).^2); %get euclidian distance to other points
            [~,b] = sort(cDist, 'ascend'); %sort to get the closest neighbours
            cID(:, iComps) = othGroups(b(1:nrNeighbors)); %get group ID for closes members
        end
        
        % check result for each component
        groupCompCnt = NaN(nrGroups, testPCs);
        for x = 1 : nrGroups
            groupCompCnt(x, :) = sum(cID == allGroups(x),1);
        end

        % classifier performance based on individual components
        [~,b] = max(groupCompCnt);
        estCompClass{iMice}(iRuns,:) = b == cGroup;
    end
    
    % compute classifier accuracy per recording.
    % this calculates the average probability of recognizing the cell type of 
    % a given recording by using a single component of this UNSEEN recording
    cRecs = unique(compRecIdx{iMice});
    tempPerf = NaN(1, length(cRecs));
    for iRecs = 1 : length(cRecs)
        tempPerf(iRecs) = mean(mean(estCompClass{iMice}(:, compRecIdx{iMice} == iRecs), 1));
    end
    decoderPerf{allGroups == cGroup} = [decoderPerf{allGroups == cGroup}, tempPerf];
    fprintf('Done %d/%d mice\n', iMice, length(cMice))
    
end


%% show classifier result
figure('renderer','painters');
plotIdx = [1 3 2];
for x = 1:length(plotIdx)
    Violin(decoderPerf{plotIdx(x)}, x, 'ViolinColor', groupColor{plotIdx(x)}, 'bandwidth', 0.05);
    [h(x), p(x)] = ttest(decoderPerf{plotIdx(x)},0.33);
    n(x) = length(decoderPerf{plotIdx(x)});
end
title('single component classifier accuracy');
ax = gca; ax.XTick = 1 : nrGroups; ax.XTickLabel = groups(plotIdx);
ylabel('percent correct');
axis square; ylim([0 1.05]); xlim([0.5 4.5]);  hline(1/3);

disp(['whole frame component classifier p-values: ' num2str(p)]);
disp(['whole frame component session count: ' num2str(n)]);


%% check for area sizes and predictive regions
a = cat(2, estCompClass{:});
a = mean(a,1);
strongAreas = a > 0.99;
weakAreas = ~strongAreas;

strongSize = sum(X(:, strongAreas)>0.2, 1);
weakSize = sum(X(:, weakAreas)>0.2, 1);
        
strongSize = sqrt(strongSize') / 51.5;
weakSize = sqrt(weakSize') / 51.5;

figure('renderer','painters');
violinPlotAlt(strongSize, 'histOri', 'left', 'widthDiv', [2 1], 'showMM', 0, 'color', mat2cell([1, 0, 0], 1));
violinPlotAlt(weakSize, 'histOri', 'right', 'widthDiv', [2 2], 'showMM', 0, 'color',  mat2cell([0, 0, 1], 1));
set(gca, 'xtick', [0.6 1.4], 'xticklabel', {'predictive', 'mixed'}, 'xlim', [0.2 1.8]);
ylabel('area size (mm^2)'); axis square
[pval, ~] = ranksum(strongSize, weakSize);
ax = gca;
mysigstar(ax, xticks, ax.YLim(2), pval);
ax.YLim(2) = ax.YLim(2) + 0.1;
ax.YTick = 0 : 0.25 : 1.5;
title(sprintf('pVal:%d; nrSp:%d; nrUSp:%d', pval, length(strongSize), length(weakSize)));
fprintf('median unspecific size: %f\n', nanmedian(weakSize))
fprintf('median preditive size: %f\n', nanmedian(strongSize))


%% make decoder performance plot for all regions
splitRegionMap = regionMap;
a = cat(2, estCompClass{:});
a = mean(a,1);

load([umapPath compFile.name], 'groupIdx', 'mouseIdx', 'areaIdx', 'recIdx', 'shrinkMask')
load([umapPath umapFile.name], 'umapOut', 'hdbLabels')
umapOut(:,1) = umapOut(:,1) - min(umapOut(:,1));
umapOut(:,2) = umapOut(:,2) - min(umapOut(:,2));
umapOut = round(umapOut * 100);

areaIdx = abs(areaIdx(~rmvIdx));
umapOut = umapOut(~rmvIdx,:);
mouseIdx = mouseIdx(~rmvIdx);
groupIdx = groupIdx(~rmvIdx);
allAreas = unique(areaIdx);
allAreas(abs(allAreas) == 1) = [];

figure('renderer','painters');
areaPerf = NaN(nrGroups, length(allAreas));
perfMap = NaN(numel(splitRegionMap), nrGroups);
compCntMap = NaN(numel(splitRegionMap), nrGroups);
for iGroups = 1 : nrGroups
    Cnt = 0;
    for iAreas = allAreas
        
        Cnt = Cnt + 1;
        cIdx = areaIdx == iAreas & groupIdx == iGroups;
        areaPerf(iGroups, Cnt) = nanmean(a(cIdx));
        perfMap(splitRegionMap(:) == iAreas, iGroups) = areaPerf(iGroups, Cnt);
        compCntMap(splitRegionMap(:) == iAreas, iGroups) = sum(cIdx);
        
    end
    
    subplot(1, nrGroups, iGroups);
    imageScale(reshape(perfMap(:,iGroups), size(allenMask)));
    colormap inferno; rateDisc_plotAllenOutline(gca);
    caxis([0 0.8]); title(groups(allGroups(iGroups)));
    
end
perfMap = reshape(perfMap, [size(allenMask), nrGroups]);



%% results for temporal components
targRegion = -2;
cFile = strrep(fullfile(compFile.folder, compFile.name), '.mat', ['_umap_area_' num2str(targRegion) '.mat']);
load([umapPath compFile.name], 'areaIdx')
areaIdx = areaIdx(~rmvIdx);

load(cFile, 'tempUmapOut', 'umapOut', 'localGroupIdx');
tempUmapOut = tempUmapOut(localGroupIdx ~= rmvGroup, :);
umapOut = umapOut(localGroupIdx ~= rmvGroup, :);

tempUmapOut(:,1) = tempUmapOut(:,1) - min(tempUmapOut(:,1));
tempUmapOut(:,2) = tempUmapOut(:,2) - min(tempUmapOut(:,2));
tempUmapOut = round(tempUmapOut * 100);

h1 = figure('renderer','painters'); hold on;
cIdx = areaIdx == targRegion;
cMice = unique(mouseIdx(groupIdx ~= rmvGroup));
for x = 1 : length(cMice)
    
    useIdx = cIdx  & mouseIdx == cMice(x);
    cGroup = unique(groupIdx(useIdx));
    arrayIdx = mouseIdx(cIdx) == cMice(x); %find current mouse in local data array
    
    scatter(tempUmapOut(arrayIdx,1), tempUmapOut(arrayIdx,2), [], groupIdx(useIdx), ...
    'Marker', Markers{x}, 'SizeData', 50, 'MarkerEdgeColor', groupColor{cGroup}); axis square; hold on;

end
ax = gca;
ax.XTick = []; ax.YTick = [];
title(['umap, all PCs: region' num2str(targRegion)]);

% show location of example PCs
segIdx = [1 0.75 1.25 0.5 1];
segFrames = cumsum(floor(segIdx * 15)); %max nr of frames per segment

clear coords
coords{1} = [756,1136]; %EMX
coords{2} = [100,1140]; %FEZ
coords{3} = [560,120]; %Plexin
nrTraces = 10;

for x = 1 : length(coords)

    figure; subplot(1,2,1);

    [a, idx] = sort(sum(abs(tempUmapOut - coords{x}), 2), 'ascend');

    localAreaIdx = find(cIdx);
    cTrace = Y(:,localAreaIdx(idx(1:nrTraces)));
    cTrace = reshape(cTrace, 66, [], nrTraces);
    % cTrace = squeeze(nanmean(cTrace(:,[1 3],:),2));
    cTrace = squeeze(nanmean(cTrace,2));
    traceStd = nanstd(cTrace, [], 1);
    cTrace = cTrace ./ traceStd;
    
    % cTrace = (cTrace(1:66,:) + cTrace(67:end,:)) / 2;
    % cTrace = (cTrace(1:66,:));
    arrayPlot(cTrace, 'color', groupColor{x}); title([groups{x} '; region' num2str(targRegion)]); niceFigure; axis square
    nvline(segFrames, 'k--'); ylabel('SDUs');

    subplot(1,2,2);
    cMaps = X(:,localAreaIdx(idx(1)));
    cMap = arrayShrink(nanmean(cMaps,2), shrinkMask);
    cImg = imageScale(cMap); caxis([0 1]);
    colormap inferno; rateDisc_plotAllenOutline(gca);

    figure(h1);
    plot(tempUmapOut(idx(1:nrTraces),1), tempUmapOut(idx(1:nrTraces),2), 'ko', 'MarkerSize', 10)

end


