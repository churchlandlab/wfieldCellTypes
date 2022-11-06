% code to analyze locaNMF results and check changes over time
% cPath = '\\churchlandNAS\homes\DOMAIN=CSHL\smusall\BpodImager\Animals\'; %path to churchlandNAS
cPath = 'X:\RawData\'; %local path
animals = {'CSP22' 'CSP23' 'CSP32' 'CSP38' ...
           'Fez71' 'Fez72' 'Fez73' 'Fez74' 'Fez75' ...
           'Plex60' 'Plex61' 'Plex65' 'Plex66' ...
           'mSM63' 'mSM64' 'mSM65' 'mSM66'};
       
load('allenDorsalMapSM.mat')
allenMask = dorsalMaps.allenMask;
trainingRange = 'allAudio'; %use this to run analysis only in a certain range of training
umapPath = 'Q:\BpodImager\umapClust\'; %save path for umap
savePath = 'Q:\BpodImager\spatialClust\'; %save path for clustering
cTime = datestr(now,'yyyymmddTHHMM');
shrinkFact = 2;
groups = {'mSM', 'Fez', 'Plex' 'CSP'};
recIdx = [];
mouseIdx = [];
groupIdx = [];
areaIdx = [];

opts.preStim = 3;
opts.frameRate = 15;
segIdx = [1 0.75 1.25 0.5 1];
segFrames = cumsum(floor(segIdx * opts.frameRate)); %max nr of frames per segment

%% get regions
[regionMap, regionLabels] = rateDisc_makeRegionMap; %get regionmap and labels

[areaMask, areaLabels] = rateDisc_areaMasks(allenMask); %get different allen areas
shrinkMask = allenMask;
% cMask = areaMask{strcmpi(areaLabels, 'MOB')};  %don't use OB areas
% shrinkMask(cMask) = true;
shrinkMask = ~bwareaopen(~imclose(shrinkMask,strel('disk',4)), 100);
shrinkMask = arrayResize(shrinkMask,shrinkFact) == 1;
shrinkMask(:, round(size(shrinkMask,2)/2)-2 : round(size(shrinkMask,2)/2)+2) = true; %remove midline

%% run over animals
allA = cell(1,length(animals)); %spatial components
allC = cell(1,length(animals)); %spatial components
trialCnt = cell(1,length(animals)); %keep trialcount for each recording
aCnt = 0;
for iAnimals = 1 : length(animals)
    aCnt = aCnt + 1;
    rCnt = 0;
    
    %current animal
    cAnimal = animals{iAnimals}; % current animal
    load([cPath cAnimal filesep 'SpatialDisc' filesep 'recs_allAudio'])

    fprintf('Current animal: %s\n', cAnimal);
    
    %go through recordings
    allA{iAnimals} = cell(1, length(recs));
    allC{iAnimals} = cell(1, length(recs));
    trialCnt{iAnimals} = NaN(1,length(recs));
    for iRecs = 1 : size(recs,1)
        
        cRec = strtrim(regexprep(recs(iRecs,:),char(0),''));
        fPath = [cPath cAnimal filesep 'SpatialDisc' filesep cRec filesep]; %Widefield data path
        
        try
            load([fPath 'newAC_20_50.mat'], 'A', 'C', 'areas');
            A = A(:,:,areas ~= 1 & areas ~= -1); %don't use OB areas
            C = C(areas ~= 1 & areas ~= -1,:); %don't use OB areas
            areas = areas(areas ~= 1 & areas ~= -1); %don't use OB areas
            A = arrayResize(A, shrinkFact);
            A = arrayCrop(A, shrinkMask);
            A = reshape(A, [], size(A,3));
            for x = 1 : size(A,2)
                A(isnan(A(:,x)),x) = prctile(A(:,x), 10);
            end
            
            A = reshape(A, size(shrinkMask,1), size(shrinkMask,2), []);
            for x = 1 : size(A,3)
                cImg = A(:,:,x);
                cImg = imadjust(cImg,[min(cImg(:)) prctile(cImg(:), 99.9)], [], 2);
                cImg = arrayFilter(cImg, 2, 2, 2);
                cImg = mat2gray(cImg);
                A(:,:,x) = cImg;
            end
            
            % keep spatial components and trialcount
            allA{iAnimals}{iRecs} = arrayShrink(A, shrinkMask, 'merge');
            load([fPath 'Vc.mat'], 'bTrials');
            trialCnt{iAnimals}(iRecs) = length(bTrials);
            areaIdx = [areaIdx, areas];
            
            % align temporal components to task episodes
            try
                load([fPath 'rsVc.mat'], 'Vc', 'bTrials');
                if size(Vc,2) * size(Vc,3) < size(C,2)
                    error;
                end
            catch
                load([fPath 'Vc.mat'], 'Vc', 'bTrials');
            end
            bhvFile = dir([fPath cAnimal '_SpatialDisc*.mat']);
            load([fPath bhvFile(1).name], 'SessionData');
                
            [~,A,B] = size(Vc);
            Vc = reshape(Vc, size(Vc,1), []);
            newC = NaN(size(C,1), size(Vc,2), 'single');
            cIdx = ~isnan(Vc(1,:));
            newC(:, cIdx) = C;
            newC = reshape(newC, size(C,1), A, B);
            bhv = selectBehaviorTrials(SessionData,bTrials); %only use completed trials that are in the Vc dataset
            newC = rateDisc_getBhvRealignment(newC, bhv, segFrames, opts); %align to different trial segments
            
            %split choice sides for a bit more variance
            temp{1} = nanmean(newC(:,:,bhv.ResponseSide == 1),3);
            temp{2} = nanmean(newC(:,:,bhv.ResponseSide == 2),3);
            temp{3} = nanmean(newC(:,:,bhv.CorrectSide == 1),3);
            temp{4} = nanmean(newC(:,:,bhv.CorrectSide == 2),3);
            allC{iAnimals}{iRecs} = cat(2,temp{:}); %keep average over all aligned trials for each component

            fprintf('Recording %d/%d\n', iRecs, length(recs));
            clear A bTrials
            
            %make indicies
            rCnt = rCnt + 1;
            cRec = ones(1,size(allA{iAnimals}{iRecs},2), 'single') .* rCnt;
            recIdx = [recIdx, cRec];
            
            cRec = ones(1,size(allA{iAnimals}{iRecs},2), 'single') .* aCnt;
            mouseIdx = [mouseIdx, cRec];
            
            cGroup = find(contains(groups, animals{iAnimals}(1:3)));
            cGroup = ones(1,size(allA{iAnimals}{iRecs},2), 'single') .* cGroup;
            groupIdx = [groupIdx, cGroup];
            
        catch ME
            disp(ME.message);
        end
    end
end

%% merge clusters into one large array.
% spatial clusters
X = cat(2,allA{:});
X = single(cat(2, X{:}));

% temporal clusters
Y = cat(2,allC{:});
Y = single(cat(1, Y{:}));

% save data
save([umapPath cTime '_localA.mat'], 'X', 'Y', 'shrinkMask', 'recIdx', 'mouseIdx', 'groupIdx', 'areaIdx', 'trialCnt', 'animals', '-v7.3');
