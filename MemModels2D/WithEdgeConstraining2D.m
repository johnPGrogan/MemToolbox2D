% WithEdgeConstraining2D(model)
% This wrapper function can be used to adjust the bivariate distributions
% on items such that responses falling outside the screen are shifted to
% the nearest point on the edge. We integrate the distribution that falls
% outside the screen and add that to the responses falling on the edge.
% This mirrors the generative process when boundary=1.
% 
% You need to have the 'targets' field in your data structure:
%   data.targets : the raw target locations [x;y] on each trial (for this
%                   condition)
% 
% This 2D version is adapted from code in the MemToolbox package(Suchow, J. W.,
% Brady, T. F., Fougnie, D., & Alvarez, G. A. (2013). Modeling visual working 
% memory with the MemToolbox. Journal of Vision, 13(10):9, 1–8. 
% doi:10.1167/13.10.9. MemToolbox.org). 
% It was adapted by John P. Grogan, 2019.
%

function model = WithEdgeConstraining2D(model, constrain)

% check that the model is a 2D model
if ~regexp('2D',model.name)
    warning('you are using WithEdgeConstraining2D() on a model which does not seem to be a 2D model')
end

if ~exist('constrain','var') % constrain by default
	constrain = 1;
end

% Take model and turn it into a model with a bias term
model.name = [model.name ' with edge constraining in pdf'];

% Adjust pdf and prior
model.oldPdf = model.pdf;
model.pdf = @NewPDF;

    function [p, pSep] = NewPDF(data, varargin)
    % adjust likelihoods by proportion of distributions falling within screen
        
        [~,pSep] = model.oldPdf(data, varargin{:}); % get separate likelihoods for error types
        
        if constrain
            if ~isfield(data, 'targets')
                error('data structure must contain targets coordinates');
            end
            
            % calculate distance from target to far edge (other edge is [0;0], so -data.targets)
            if ~isfield(data, 'edgeDists')
                data.edgeDists = data.dimensions - data.targets;
            end
            
            % renormalise pdf by the proportion of responses falling inside
            % screen
            ratioOutsideScreen = mvnpdfRatio(data, varargin{end});
            
            % find which columns to update (not guessing)
            guessInd = ~cellfun(@isempty,regexp('g', model.paramNames));
            
            pSep(:,~guessInd) = pSep(:,~guessInd)  + ratioOutsideScreen; % normalise
            
        end

        p = sum(pSep,2); % sum across diff response types to get p per trial
    
    end


    function ratioOutside= mvnpdfRatio(data, sigma)
        % function ratioOutside = mvnpdfRatio(data, sigma)
        % Calculate the proportion of a mvnpdf that falls outside the
        % screen, within a certain area. This is added to the pdf for
        % responses falling on screen edges, assuming that people are
        % simply shifting any generated responses falling outside to lie on
        % the nearest edge. 
        % 
        % Inputs: 
        %       data = MemToolbox2D data structure, must contain data.targets, and
        %               data.dimensions
        %       sigma = sigma (SD) value for mvnpdf
        % 
        % Outputs: 
        %       propOutside = ratio of integral outside screen where the
        %       response is the closest valid point, to the entire integral


        nTrials = size(data.errors,2); 
        
        % copied from mvnpdf.m, made faster by removing excess stuff
        mvnFunc = @(x,y) reshape(exp(-0.5*sum ([col(x),col(y)].^2 ./ [sigma sigma].^2,2) - sum(log(sqrt([sigma sigma].^2)),2) - 2*log(2*pi)/2), size(x));
        % slower version
        % mvnFunc = @(x,y) reshape(mvnpdf([col(x),col(y)], [0,0], [sigma^2,sigma^2]), size(x));


        % get whether each response is on edge or not
        % -1=left/bottom, +1=right/top
        edges = [(data.errors(1,:) <= -data.targets(1,:))*-1 + (data.errors(1,:) >= data.edgeDists(1,:));
                    (data.errors(2,:) <= -data.targets(2,:))*-1 + (data.errors(2,:) >= data.edgeDists(2,:));];
        

        % boxes to integrate over outside the screen
        allBoxes = data.errors([1 1 2 2],:)'; % preset these
        allBoxes(all(edges==0,1)',:) = 0; % set insides to nan
        
        width = 1; % with of strip to integrate outside screen
        widths = zeros(size(allBoxes,1),4); % preset
        
        
        % if on edge, change limits
        allBoxes(edges(1,:)==-1, 1) = -data.dimensions(1); % left
        allBoxes(edges(1,:)==1, 2) = data.dimensions(1); % right
        
        widths(edges(1,:)~=0,3:4) = repmat([-width/2, width/2],sum(edges(1,:)~=0),1); % if left/right, set width on y axis
        
        % now top/bottom
        allBoxes(edges(2,:)==-1, 3) = -data.dimensions(2); % bottom
        allBoxes(edges(2,:)==1, 4) = data.dimensions(2); % top
        
        widths(edges(2,:)~=0,1:2) = repmat([-width/2, width/2],sum(edges(2,:)~=0),1); % if top/bottom set width on x axis
        
        % if corner, no width
        widths(all(edges~=0,1),:) = 0;
        
        allBoxes = allBoxes + widths; % add the width
        
        isOnEdge = find(any(edges~=0,1)); % find trials to adjust

        % integrate each of these boxes
        AOut = zeros(1, nTrials);
        for j = 1:length(isOnEdge) % for each outside point
            i = isOnEdge(j);
            AOut(1,i) = integral2(mvnFunc, allBoxes(i,1),allBoxes(i,2),allBoxes(i,3),allBoxes(i,4));
        end

        % over max area allowed
        ATot = integral2(mvnFunc, -data.dimensions(1), data.dimensions(1), -data.dimensions(2), data.dimensions(2));
        
        ratioOutside = (AOut ./ ATot)'; % divide pdf by this to adjust for portion of pdf outside of screen
    end

end
