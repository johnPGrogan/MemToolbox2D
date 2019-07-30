% GENERATEDISPLAYS2D can be called to generate 2D locations to present in a working memory task.
%
%  displays = GenerateDisplays2D(numTrials, itemsPerTrial, mode, dimensions, minDists)
%
% numTrials must be an integer. itemsPerTrial can be either an integer or
% a vector the 1 x numTrials in size, specifying the number of items on
% each trial separately.
%
% 'mode' chooses the method by which the colors are selected.
% Currently the only supported mode is drawing randomly from a multivariate
% normal distribution.
%    mode = 1: returns information needed to do continuous report
% 	 mode = 2: returns info for 2AFC 
% 
% dimensions is optional which sets the screen dimensions in whatever units
% are being used. Defaults is 1366*768 (i.e. pixels)
% 
% minDists is a vector of minimum distances between items in the same trial
% and between items and edges and between items and screen centre. 
% If only one value is supplied, it is used for all distances. Default is 0
%
% This 2D version is adapted from code in the MemToolbox package(Suchow, J. W.,
% Brady, T. F., Fougnie, D., & Alvarez, G. A. (2013). Modeling visual working 
% memory with the MemToolbox. Journal of Vision, 13(10):9, 1–8. 
% doi:10.1167/13.10.9. MemToolbox.org). 
% It was adapted by John P. Grogan, 2019.
%

function displays = GenerateDisplays2D(numTrials, itemsPerTrial, mode, dimensions, minDists)

% Default mode is 1 for continuous error
if nargin < 3
    mode = 1;
end
if nargin < 4
    dimensions = [1366; 768];
end
if nargin < 5
    minDists = [0 0 0];
elseif length(minDists)==1
    minDists = repmat(minDists,1,3);
elseif length(minDists) ~= 3
    error('minDists should be length 1 or 3')
end
    
if mode == 1 || mode == 2
    % Generate random items - x & y coords for each item
    for i = 1:max(itemsPerTrial)
        displays.items((2*i-1):2*i,:) = [unifrnd(0, dimensions(1), 1, numTrials);unifrnd(0, dimensions(2), 1, numTrials)];
        displays.items((2*i-1):2*i,i>itemsPerTrial) = NaN;%remove unused items on each trial
    end
    
    % if minDists are used, replace any trials that are within that
    % distance
    if any(minDists > 0) % if there is a minimum item distance
        toReplace = GetItemDists(displays.items, dimensions, minDists); % which trials have items to close to others or edges
        nRepIter = 1;
        while any(toReplace)
            for i = 1:max(itemsPerTrial)
                displays.items((2*i-1):2*i,toReplace) = [unifrnd(0, dimensions(1), 1, sum(toReplace));unifrnd(0, dimensions(2), 1, sum(toReplace))];
                displays.items((2*i-1):2*i,i>itemsPerTrial) = NaN;%remove unused items on each trial
            end        
            toReplace = GetItemDists(displays.items, dimensions, minDists); % get distances between each item on each trial
            nRepIter = nRepIter + 1;
        end
    end
    
    displays.whichIsTestItem = ceil(rand(1, numTrials).*itemsPerTrial); % select chosen item
    
    for i=1:size(displays.items,2) % get target coords
        whichTest = displays.whichIsTestItem(i);%get which item was the target
        displays.targets(:,i) = displays.items((whichTest*2-1):whichTest*2,i);
    end
    
    displays.dimensions = dimensions;
    % Extract useful information for models
    displays = AddUsefulInfo(displays);
    
    % If 2AFC, choose a change size for each tested item also
    if mode==2
        % pick random coords of moved target
        displays.changeSize(1,:) = unifrnd(0, dimensions(1), 1, numTrials);
        displays.changeSize(2,:) = unifrnd(0, dimensions(2), 1, numTrials);
        
        displays.changeSize = displays.changeSize - displays.targets; % get distance from targets
    end
    
else
    warning('No such mode, defaulting to 2D distance.')
    displays = GenerateDisplays2D(numTrials,itemsPerTrial,1,dimensions);
end
end

function displays = AddUsefulInfo(displays)
% Add in the distance of the distractors to the target for each trial
    
    allItems = 1:size(displays.items,1)/2; % number of items, 2 rows per item
    
    nDist = max(allItems) - 1;%number of distractors
    
    for i=1:size(displays.items,2) % for each trial
        
        whichTest = displays.whichIsTestItem(i);%get which item was the target

        distInds = allItems(allItems~=whichTest);%get indices of row-pairs for distractors
        distInds = sort([distInds*2-1, distInds*2]);%sort
        
        %get distance of target to all distractors for this trial
        displays.distractors(:,i) = displays.items(distInds, i) - ...
            repmat(displays.items((whichTest*2-1):whichTest*2,i),nDist,1);
    end
    
    % Add set size for each trial
    displays.n = sum(~isnan(displays.items)) ./ 2;
   
end

function toReplace = GetItemDists(items,dimensions, minDists)
% toReplace = GetItemDists(items,dimensions, minDists)
% get euclidean distances between all items on each trial, and distance
% from edges and centre. return logical vector of which trials have items
% within those minDists
% items is displays.items

    [nRows,nCols] = size(items);

    itemComplex = complex(items(1:2:nRows,:), items(2:2:nRows,:)); % make into complex coords ox x(real) and y (imag)
    nItems = size(itemComplex,1); % number of items

    itemDists = NaN(nchoosek(nItems,2),nCols); % preallocate for all pairs
    k = 1;
    for i = 1:nItems
        for j = i+1:nItems
            if i~=j % for diff items
                itemDists(k,:) = abs(itemComplex(i,:) - itemComplex(j,:));
                k = k + 1;
            end
        end
    end

    edgeDists = NaN(nItems, nCols);
    edgeDists(1,:) = min(real(itemComplex),[],1); % x coord is distance from x=0
    edgeDists(2,:) = min(imag(itemComplex),[],1); % same for y coord
    %distance from upper bounds
    edgeDists(3:4,:) = dimensions - [max(real(itemComplex),[],1);  max(imag(itemComplex),[],1)];
    
    centreDists = NaN(nItems, nCols);
    centre = complex(dimensions(1)./2, dimensions(2)./2); % complex coords of screen centre
    for i = 1:nItems
        centreDists(i,:) = abs( itemComplex(i,:) - centre ); %distance from centre
    end
    toReplace = min(itemDists,[],1) < minDists(1) | ...
                min(edgeDists,[],1) < minDists(2) | ...
                min(centreDists,[],1) < minDists(3);


end
