% code to analyze umap clustering whole frame locaNMF components
load('allenDorsalMapSM.mat')
allenMask = dorsalMaps.allenMask;
groups = {'mSM' 'Fez' 'Plex' 'CSP'};
cPath = 'X:\RawData\'; %local path
clustPath = 'G:\spatialClust\';
nrNeighbors = 10; %number of neighbors for classifier
nrRuns = 100; %nr of runs when testing each individual component
pcCnt = 500; %maximum nr of components from each animal per run
noCStr = true;
groupColor = {[1 0 0] [0 0 1] [0 1 0] [0 1 1]}; %colors for groups
rng(0);

umapPath = 'X:\RawData\Celltypes_dataset\umap\';
umapFilename = '20210309T2356_fullA';

%% show locanMF regions
load([cPath 'localRegionMap.mat'])
regionMap = localRegionMap;
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
colormap(cImg.Parent,hsv(256));
caxis([0 max(regionMapPlot(:))]); hold on;
for iRegions = (1:10) + max(regionMap(:))
    cFrame = regionMapPlot == iRegions;
    areaOutline= bwboundaries(cFrame); %outline of selected area
    for x = 1 : length(areaOutline)
        plot(areaOutline{x}(:,2),areaOutline{x}(:,1),'k')
    end
end


%% show umap result
% get all clusters
if noCStr
    rmvGroup = 4; %dont use CStr neurons
else
    rmvGroup = 0; %use CStr neurons
end

umapFile = dir([umapPath '*' umapFilename '_umap.mat']); % load this for wholeFrame - sNMF components
compFile = dir([umapPath '*' umapFilename '.mat']); % load this for wholeFrame - sNMF components

load([umapPath umapFile.name], 'umapOut', 'hdbLabels')
load([umapPath compFile.name], 'groupIdx', 'mouseIdx', 'X', 'shrinkMask')

rmvIdx = groupIdx == rmvGroup;
umapOut(:,1) = umapOut(:,1) - min(umapOut(:,1));
umapOut(:,2) = umapOut(:,2) - min(umapOut(:,2));
umapOut = round(umapOut * 100);

umapOut = umapOut(~rmvIdx,:);
X = X(:,~rmvIdx);
mouseIdx = mouseIdx(~rmvIdx);
groupIdx = groupIdx(~rmvIdx);


%% show UMAP clusters for whole frame components
h1 = figure('renderer','painters'); hold on;
Markers = {'+','o','*','x','v','d','^','s','>','<','+','o','*','x','v','d','^','s','>','<'};
cMice = unique(mouseIdx);
cMice = cMice(randperm(length(cMice)));
for x = 1 : length(cMice)
    
    useIdx = mouseIdx == cMice(x);
    cGroup = unique(groupIdx(useIdx));
    
    scatter(umapOut(useIdx,1), umapOut(useIdx,2), [], groupIdx(useIdx), ...
    'Marker', Markers{x}, 'SizeData', 20, 'MarkerEdgeColor', groupColor{cGroup}); axis square;
end
ax = gca;
ax.XTick = []; ax.YTick = [];
title('umap, all PCs');


%% plot EMX subset with individal mice - first PCs
yIdx = umapOut(:,2) > 610 & umapOut(:,2) < 675;
xIdx = umapOut(:,1) > 0 & umapOut(:,1) < 150;
cIdx = xIdx & yIdx;

figure(h1);
cMice = unique(mouseIdx(cIdx));
cMice = cMice(randperm(length(cMice)));
for x = 1 : length(cMice)
    
    useIdx = cIdx'  & mouseIdx == cMice(x);
    cGroup = unique(groupIdx(useIdx));
    
    scatter(umapOut(useIdx,1), umapOut(useIdx,2), [], groupIdx(useIdx), ...
    'Marker', Markers{x}, 'SizeData', 100, 'MarkerEdgeColor', groupColor{cGroup}); axis square;
end
ax = gca;
ax.XTick = []; ax.YTick = [];
title('umap, all PCs');

% show location of example PCs
h2 = figure('renderer','painters'); hold on;
selectIdx = [5 10 15 45];
fezPCs = find(cIdx);
for x = 1 : 4

    figure(h2);
    subplot(2,2,x);
    cMap = X(:,fezPCs(selectIdx(x)));
    cMap = arrayShrink(cMap, shrinkMask);
    imageScale(cMap); axis image; colormap inferno; drawnow;
    caxis([0 prctile(cMap(:), 95)]);    

    figure(h1); hold on;
    plot(umapOut(fezPCs(selectIdx(x)),1), umapOut(fezPCs(selectIdx(x)),2), 'ko', 'LineWidth', 3);    
end

%% plot Plexin subset with individal mice - first PCs
yIdx = umapOut(:,2) > 520 & umapOut(:,2) < 620;
xIdx = umapOut(:,1) > 150 & umapOut(:,1) < 340;
cIdx = xIdx & yIdx;

figure(h1);
cMice = unique(mouseIdx(cIdx));
for x = 1 : length(cMice)
    
    useIdx = cIdx'  & mouseIdx == cMice(x);
    cGroup = unique(groupIdx(useIdx));
    
    scatter(umapOut(useIdx,1), umapOut(useIdx,2), [], groupIdx(useIdx), ...
    'Marker', Markers{x}, 'SizeData', 100, 'MarkerEdgeColor', groupColor{cGroup}); axis square;
end
ax = gca;
ax.XTick = []; ax.YTick = [];
title('umap, all PCs');

% show location of example PCs
h2 = figure('renderer','painters'); hold on;
selectIdx = [20 50 85 122];
fezPCs = find(cIdx);
for x = 1 : 4

    figure(h2);
    subplot(2,2,x);
    cMap = X(:,fezPCs(selectIdx(x)));
    cMap = arrayShrink(cMap, shrinkMask);
    imageScale(cMap); axis image; colormap inferno; drawnow;
    caxis([0 prctile(cMap(:), 95)]);    

    figure(h1); hold on;
    plot(umapOut(fezPCs(selectIdx(x)),1), umapOut(fezPCs(selectIdx(x)),2), 'ko', 'LineWidth', 3);    
end

%% plot Fez subset with individal mice - first PCs
yIdx = umapOut(:,2) > 300 & umapOut(:,2) < 400;
xIdx = umapOut(:,1) > 350 & umapOut(:,1) < 490;
cIdx = xIdx & yIdx;

figure(h1);
cMice = unique(mouseIdx(cIdx));
for x = 1 : length(cMice)
    
    useIdx = cIdx'  & mouseIdx == cMice(x);
    cGroup = unique(groupIdx(useIdx));
    
    scatter(umapOut(useIdx,1), umapOut(useIdx,2), [], groupIdx(useIdx), ...
    'Marker', Markers{x}, 'SizeData', 100, 'MarkerEdgeColor', groupColor{cGroup}); axis square;
end
ax = gca;
ax.XTick = []; ax.YTick = [];
title('umap, all PCs');

% make some example PCs
h2 = figure('renderer','painters'); hold on;

selectIdx = [1 25 50 75];
fezPCs = find(cIdx);
for x = 1 : 4

    figure(h2);
    subplot(2,2,x);
    cMap = X(:,fezPCs(selectIdx(x)));
    cMap = arrayShrink(cMap, shrinkMask);
    imageScale(cMap); axis image; colormap inferno; drawnow;
    caxis([0 prctile(cMap(:), 95)]);    

    figure(h1); hold on;
    plot(umapOut(fezPCs(selectIdx(x)),1), umapOut(fezPCs(selectIdx(x)),2), 'ko', 'LineWidth', 3);    
end

%% plot Fez subset with individal mice - other PCs
yIdx = umapOut(:,2) > 625 & umapOut(:,2) < 670;
xIdx = umapOut(:,1) > 1375 & umapOut(:,1) < 1400;
cIdx = xIdx & yIdx;

figure(h1);
cMice = unique(mouseIdx(cIdx));
for x = 1 : length(cMice)
    
    useIdx = cIdx'  & mouseIdx == cMice(x);
    cGroup = unique(groupIdx(useIdx));
    
    scatter(umapOut(useIdx,1), umapOut(useIdx,2), [], groupIdx(useIdx), ...
    'Marker', Markers{x}, 'SizeData', 100, 'MarkerEdgeColor', groupColor{cGroup}); axis square;
end
ax = gca;
ax.XTick = []; ax.YTick = [];
title('umap, all PCs');

% make some example PCs
h2 = figure('renderer','painters'); hold on;

selectIdx = [4 15 18 32];
fezPCs = find(cIdx);
for x = 1 : 4

    figure(h2);
    subplot(2,2,x);
    cMap = X(:,fezPCs(selectIdx(x)));
    cMap = arrayShrink(cMap, shrinkMask);
    imageScale(cMap); axis image; colormap inferno; drawnow;
    caxis([0 prctile(cMap(:), 95)]);    

    figure(h1); hold on;
    plot(umapOut(fezPCs(selectIdx(x)),1), umapOut(fezPCs(selectIdx(x)),2), 'ko', 'LineWidth', 3);    
end


%% new cell type classifier - based on umap projections from other mice
baseName = [umapFilename '_mouseTest']; %base name for umap files
load([umapPath [baseName '1.0.mat']], 'mouseIdx', 'groupIdx')

allGroups = unique(groupIdx);
cMice = unique(mouseIdx(groupIdx ~= rmvGroup));
estCompClass = cell(1,length(cMice));
compRecIdx = cell(1,length(cMice));
decoderPerf = cell(1, length(allGroups));
for iMice = 1 : length(cMice)
    
    cMouse = cMice(iMice);
    cFile = [baseName num2str(cMouse) '.0.mat'];
    load([umapPath cFile], 'umapRed', 'umapTest', 'recIdx' , 'mouseIdx', 'groupIdx')

    compRecIdx{iMice} = recIdx(mouseIdx == cMouse); %keep IDs for each moouse to make sure you can find them again

    cGroup = unique(groupIdx((mouseIdx == cMouse)));
    testPCs = size(umapTest, 1);

    groupIdx = groupIdx(mouseIdx ~= cMouse);
    useIdx = groupIdx ~= rmvGroup; %dont use rmvGroup
    allGroups = unique(groupIdx(useIdx));
    othGroupIdx = groupIdx;
    othGroupIdx(~useIdx) = NaN;
    
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
        groupCompCnt = NaN(length(allGroups), testPCs);
        for x = 1 : length(allGroups)
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
    decoderPerf{allGroups == cGroup} = [
        decoderPerf{allGroups == cGroup}, tempPerf];
end

%% show classifier results
figure('renderer','painters');
plotIdx = [1 3 2]; clear h p n;
for x = 1:length(plotIdx)
    Violin(decoderPerf{plotIdx(x)}, x, 'ViolinColor', groupColor{plotIdx(x)}, 'bandwidth', 0.05);
    [h(x), p(x)] = ttest(decoderPerf{plotIdx(x)},0.33);
    n(x) = length(decoderPerf{plotIdx(x)});
end
title('whole frame component classifier accuracy');
ax = gca; ax.XTick = 1 : length(allGroups); ax.XTickLabel = groups(plotIdx);
ylabel('percent correct');
axis square; ylim([0 1.05]); xlim([0.5 4.5]);  hline(1/3);
disp(['whole frame component classifier p-values: ' num2str(p)]);
disp(['whole frame component session count: ' num2str(n)]);




%% nested functions
function [outIdx, minCnt] = rateDisc_animalIndex(cMouseIdx, maxCnt)
% function to return a logical index 'outIdx' with an even amount of entries 
% for all animals in 'cMouseIdx'. maxCnt can be used to reduce the number
% of entries per animal.

allMice = unique(cMouseIdx(~isnan(cMouseIdx)));

for iMice = 1 : length(allMice)
    mCnt(iMice) = sum(cMouseIdx == allMice(iMice));
end
minCnt = min([mCnt, maxCnt]); %minimum nr of entries that is available in all mice

outIdx = false(1, length(cMouseIdx));
for iMice = 1 : length(allMice)
    cIdx = find(cMouseIdx == allMice(iMice));
    outIdx(cIdx(randperm(length(cIdx), minCnt))) = true;
end
end