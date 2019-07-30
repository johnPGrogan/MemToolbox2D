% WithResponseSampling2D(model, bandwidth)
% This wrapper function can be used similarly to WithBias2D(), to make the
% guess distribution be sampled from a pdf rather than a uniform
% distribution as standard.
% 
% ksdensity is used with a normal distribution and default binwidth. If the
% data.resps values fall outside the screen dimensions (as may happen when
% simulating responses using GenerateDisplays2D), it will not constrain the
% smoothing to the screen dimensions.
%
% You can supply a fixed bandwidth (scalar or two-element vector), rather
% than using Matlab's default estimated bandwidth.
% 
% You need to have two extra fields in your data structure:
%   data.resps : the raw response locations [x;y] on each trial (for this
%                condition
%   data.respPdf : a cell containing the response pdf you wish to use for 
%                  the guesses, which takes the form of a sample 
%                  distribution e.g. all the responses for a person across 
%                  all conditions and trials.
% 
% This 2D version is adapted from code in the MemToolbox package(Suchow, J. W.,
% Brady, T. F., Fougnie, D., & Alvarez, G. A. (2013). Modeling visual working 
% memory with the MemToolbox. Journal of Vision, 13(10):9, 1–8. 
% doi:10.1167/13.10.9. MemToolbox.org). 
% It was adapted by John P. Grogan, 2019.
%

function model = WithResponseSampling2D(model, bandwidth)
% If no prior is specified, default to an improper uniform prior

% check that the model is a 2D model
if ~regexp('2D',model.name)
    warning('you are using WithResponseSampling2D() on a model which does not seem to be a 2D model')
end

if ~exist('bandwidth','var') % if no bandwidth is supplied, use default from Matlab's estimator
	useBandwidth = 0;
else
	useBandwidth = 1;
end

% Take model and turn it into a model with a bias term
model.name = [model.name ' with response sampling for guesses'];

% Adjust pdf and prior
model.oldPdf = model.pdf;
model.pdf = @NewPDF;

% Shift errors and/or changeSize
    function [p, pSep] = NewPDF(data, varargin)
        dimensions = GetModelDims(data); % get the model dimensions
        
        [~,pSep] = model.oldPdf(data, varargin{:}); % get separate likelihoods for error types
        if all(pSep==0,'all') % if params exceed boundaries
            p = pSep(:,1);
        else

            % we need the coordinates of responses
            if ~isfield(data,'resps')
                if isfield(data,'targets')
                    data.resps = data.errors - data.targets;
                elseif isfield(displayInfo, 'items') && isfield(displayInfo, 'whichIsTestItem')
                    whichItems = displayInfo.whichIsTestItem;
                    for i = 1:size(displayInfo.items,2)
                        targs(:,i) = displayInfo.items([whichItems(i)*2-1, whichItems(i)*2],i);
                    end
                    data.resps = data.errors - targs;
                end
            end
        
            % get response pdf prob for responses
            if any(data.resps <= 0,'all') || any(data.resps >= dimensions, 'all')
				if useBandwidth % if a set bandwidth is specified
					pG = (varargin{1}) .* ksdensity(data.respPdf{1}',data.resps','Bandwidth',bandwidth); % eval pdf for each resp
				else
					pG = (varargin{1}) .* ksdensity(data.respPdf{1}',data.resps'); % eval pdf for each resp
				end
            else
				if useBandwidth
					pG = (varargin{1}) .* ksdensity(data.respPdf{1}',data.resps','Support',[-1 -1;dimensions'],'Bandwidth',bandwidth);%eval pdf for each 
				else
					pG = (varargin{1}) .* ksdensity(data.respPdf{1}',data.resps','Support',[-1 -1;dimensions']);%eval pdf for each 
				end
            end
			guessInd = find(~cellfun(@isempty,regexp('g', model.paramNames)));
            pSep(:,guessInd+1) = pG; % replace standard guess rate

            p = sum(pSep,2); % sum for new likelihood
        end
    end

end
