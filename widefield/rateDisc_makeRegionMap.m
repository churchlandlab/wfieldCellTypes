function [regionMap, regLabels, regionAreas] = rateDisc_makeRegionMap
% function to create a map of larger cortical regions, based on the allen
% CCF. This is to simplify correlation analysis or provide a broader seed
% region for LocaNMF.

%% labels for different regions
regLabels = {'OB' 'M2a' 'M2p' 'M1' 'SSm' 'SSb' 'SSf' 'Aud' 'PPC' 'RS' 'V2' 'V1'}; 

%areas to be included in each region
regionAreas = {{'MOB'}, ...
            {'MOs', 'PL', 'FRP',}, ...
            {'MOs', 'ACAd'}, ...
            {'MOp'}, ...
            {'SSp-m'}, ...
            {'SSp-ll','SSp-ul', 'SSp-un'}, ...
            {'SSp-n', 'SSp-bfd'}, ...
            {'SSs', 'AUDd','VISal', 'VISl'}, ...
            {'SSp-tr', 'VISrl', 'VISa'}, ...
            {'RSPagl', 'RSPd', 'RSPv'}, ...
            {'VISpm', 'VISam'}, ...
            {'VISp'} ...
            };
            
%get different allen areas
[areaMask, areaLabels] = rateDisc_areaMasks('allenMask');
regionMap = NaN(size(areaMask{1}), 'single');
for x = 1 : length(regionAreas)
    
    cIdx = find(ismember(areaLabels, regionAreas{x}));
    for y = 1 : length(cIdx)
        if strcmpi(regLabels{x}, 'M2a') && strcmpi(areaLabels{cIdx(y)}, 'MOs')
            mask1 = areaMask{cIdx(y)};
            mask2 = createCircle(size(regionMap), size(regionMap)/2, 100);
            regionMap(mask1 & ~mask2) = x;

        elseif strcmpi(regLabels{x}, 'M2p') && strcmpi(areaLabels{cIdx(y)}, 'MOs')
            mask1 = areaMask{cIdx(y)};
            mask2 = createCircle(size(regionMap), size(regionMap)/2, 100);
            regionMap(mask1 & mask2) = x;            
            
        else
            regionMap(areaMask{cIdx(y)}) = x;
        end
    end
    
    % connect areas within each region
    cMask = imclose(regionMap == x, strel('disk', 8));
    regionMap(cMask) = x;
end
    
    



