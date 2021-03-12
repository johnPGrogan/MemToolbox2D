% WithEdgeResampling2D(model)
% This wrapper function can be used to normalise the bivariate
% distributions around items by the proportion of the distribution that
% falls within the the screen dimensions. e.g. if an item falls near a
% corner and only 1/4 of the pdf is valid, then that is normalised so that
% it sums to 1, thus increasing the likelihood of responses falling near
% that item.
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

function model = WithEdgeResampling2D(model, resample)

% check that the model is a 2D model
if ~regexp('2D',model.name)
    warning('you are using WithEdgeResampling2D() on a model which does not seem to be a 2D model')
end

if ~exist('resample','var') % resample by default
	resample = 1;
end

% Take model and turn it into a model with a bias term
model.name = [model.name ' with edge resampling in pdf'];

% Adjust pdf and prior
model.oldPdf = model.pdf;
model.pdf = @NewPDF;

    function [p, pSep] = NewPDF(data, varargin)
    % adjust likelihoods by proportion of distributions falling within screen
        
        [~,pSep] = model.oldPdf(data, varargin{:}); % get separate likelihoods for error types
        
        if resample
            if ~isfield(data, 'targets')
                error('data structure must contain targets coordinates');
            end
            
            % calculate distance from target to far edge (other edge is [0;0], so -data.targets)
            if ~isfield(data, 'edgeDists')
                data.edgeDists = data.dimensions - data.targets;
            end
            
            % renormalise pdf by the proportion of responses falling inside
            % screen
            ratioInsideScreen = mvnpdfRatio(data, varargin{end});
            
            % find which columns to update (not guessing)
            guessInd = ~cellfun(@isempty,regexp('g', model.paramNames));
            guessInd = [1:length(model.paramNames)] == (find(guessInd)+1);% first column is A            
            
            
            pSep(:,~guessInd) = pSep(:,~guessInd) ./ ratioInsideScreen; % normalise
            
        end

        p = sum(pSep,2); % sum across diff response types to get p per trial
    
    end


    function ratioInside = mvnpdfRatio(data, sigma)
        % function ratioInside = mvnpdfRatio(data, sigma)
        % Calculate the proportion of a mvnpdf that falls within an area
        % for the MemToolbox2D package, to adjust pdf when resampling is used.
        % 
        % Inputs: 
        %       data = MemToolbox2D data structure, must contain data.targets, and
        %               data.dimensions
        %       sigma = sigma (SD) value for mvnpdf
        % 
        % Outputs: 
        %       propInside = ratio of integral inside screen to total

%         nTrials = size(data.errors,2);
% 
%         % copied from mvnpdf.m, made faster by removing excess stuff
% %         mvnFunc = @(x,y) reshape(exp(-0.5*sum ([col(x),col(y)].^2 ./ [sigma sigma].^2,2) - sum(log(sqrt([sigma sigma].^2)),2) - 2*log(2*pi)/2), size(x));
% %         mvnFunc = @(x,y) reshape(exp(-0.5 * sum([col(x),col(y)].^2, 2) ./ sigma.^2 - 2*log(sigma) - log(2*pi)), size(x));
%         mvnFunc = @(x,y) reshape(exp(-0.5 * sum([col(x),col(y)].^2, 2) ./ sigma.^2 - log(sigma^2 * 2*pi)), size(x)); % this is fastest version, only works for [x,y]
%         % slower version
%         % mvnFunc = @(x,y) reshape(mvnpdf([col(x),col(y)], [0,0], [sigma^2,sigma^2]), size(x));
% 
%         % over max area allowed
%         ATot = integral2(mvnFunc, -data.dimensions(1), data.dimensions(1), -data.dimensions(2), data.dimensions(2));
% 
%         AIn = NaN(1, nTrials);
%         for i = 1:nTrials
% 
%             % over area within screen edges (relative to error)
%             AIn(:,i) = integral2(mvnFunc, -data.targets(1,i), data.edgeDists(1,i), -data.targets(2,i), data.edgeDists(2,i));
% 
%         end
%         
%         ratioInside = (AIn ./ ATot)'; % divide pdf by this to adjust for portion of pdf outside of screen

        dist = sum(abs((data.targets ./ data.dimensions) - .5)); % 1 = in corner, .5 = on edge, 0 = centre
        
        % get non-target locations
        nTargLoc = repmat(data.targets, size(data.distractors,1)/2, 1) + data.distractors;
        nTargDists = sum(abs((nTargLoc ./ repmat(data.dimensions, size(data.distractors,1)/2,1)) - .5)); % distance from centre
        
        f = @(x) 0.5.*x.^2 - 1.25.*x + 1; % convert from dist to prop of distribution inside (1 = all inside, 0.5 = on edge, 0.25 = corner)
        
%         plot(dist, f(dist), 'x');

%         scaling = 1 ./ f(dist); % ratio to scale by
        ratioInside = [f(dist)', f(nTargDists)'];


    end

end
