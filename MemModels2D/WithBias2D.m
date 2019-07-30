% WITHBIAS2D adds a bias terms to any 2D model
%
%  model = WithBias2D(model, dims, [priorForMu])
%
% The first parameter is the model to convert; the second (optional)
% parameter are the screen dimensions [x;y] which set the max and min
% values for the bias parameters (default is 1366*768, 
% and the third (optional) parameter is a function for the prior on mu.
% It takes two arguments (muX, muY).
%
% Thus, to take the StandardMixtureModel2D, which has a guess rate (g) and
% standard deviation (sd), and add a shift term (mu), just use:
%   model = WithBias2D(StandardMixtureModel2D())
%
% This wrapper is compatible with both and TwoAFC(). For
% example, the following works fine:
%   model = TwoAFC2D(WithBias2D(StandardMixtureModel2D));
%
% This 2D version is adapted from code in the MemToolbox package(Suchow, J. W.,
% Brady, T. F., Fougnie, D., & Alvarez, G. A. (2013). Modeling visual working 
% memory with the MemToolbox. Journal of Vision, 13(10):9, 1–8. 
% doi:10.1167/13.10.9. MemToolbox.org). 
% It was adapted by John P. Grogan, 2019.
%

function model = WithBias2D(model, dims, priorForMu )
% If no prior is specified, default to an improper uniform prior
if nargin < 3
    priorForMu = @(p)(1);
end

% check that the model is a 2D model
if ~regexp('2D',model.name)
    warning('you are using WithBias2D() on a model which does not seem to be a 2D model')
end

if ~exist('dims', 'var')
    dims = [1366; 768];
end

% Take model and turn it into a model with a bias term
model.name = [model.name ' with bias 2D'];
model.paramNames = [model.paramNames,'muX', 'muY'];
model.lowerbound = [model.lowerbound -dims(1) -dims(2)];
model.upperbound = [model.upperbound dims(1) dims(2)];
model.movestd = [model.movestd 1 1];
model.start = [model.start rand(size(model.start,1),2)*10 - 5];

% Adjust pdf and prior
model.oldPdf = model.pdf;
model.pdf = @NewPDF;

model.priorForMu = priorForMu;
if isfield(model, 'prior')
    model.oldPrior = model.prior;
    model.prior = @(p)(model.oldPrior(p(1:end-2)) .* model.priorForMu(p(end-1:end)));
end

% Adjust generator function
if isfield(model, 'generator')
    model.oldGenerator = model.generator;
    model.generator = @newGen;%(params,dims,displayInfo)(...
    %bound(model.oldGenerator(params(1:3), dims, displayInfo)+params{4} + params{5}));
end

    function y = newGen(params, dims, displayInfo)
        dimensions = GetModelDims(displayInfo); % get the model dimensions

        y = model.oldGenerator(params(1:end-2),dims,displayInfo);
        y = y + [params{end-1}; params{end}];
        
        % convert back into dimensions to apply boundary
        if isfield(displayInfo, 'targets')
            targs = displayInfo.targets;
        elseif isfield(displayInfo, 'items') && isfield(displayInfo, 'whichIsTestItem')
            whichItems = displayInfo.whichIsTestItem;
            for i = 1:size(displayInfo.items,2)
                targs(:,i) = displayInfo.items([whichItems(i)*2-1, whichItems(i)*2],i);
            end
        else
            error('need either targets or items & whichIsTestItem specified')
        end
        
        respLoc = y + targs;
        respLoc = bound(respLoc, dimensions);
        
        y = respLoc - targs;
    end

% Adjust model_plot
if isfield(model, 'modelPlot')
    model.oldModelPlot = model.modelPlot;
    model.modelPlot = @NewModelPlot;
end

    function figHand = NewModelPlot(data, params, varargin)
        dimensions = GetModelDims(data); % get the model dimensions

        if isstruct(params) && isfield(params, 'vals')
            muX = mean(params.vals(:,end-1));
            muY = mean(params.vals(:,end));
            params.vals = params.vals(:, 1:end-2);
        else
            muX = params(end-1);
            muY = params(end);
            params = params(1:end-2);
        end
        if isfield(data, 'errors')
            data.errors = bound(data.errors - [muX;muY], dimensions);
        end
        if isfield(data, 'changeSize')
            data.changeSize = bound(data.changeSize - [muX;muY],dimensions);
        end
        figHand =  model.oldModelPlot(data, params, varargin);
    end

% Shift errors and/or changeSize
    function [p, pSep] = NewPDF(data, varargin)
        dimensions = GetModelDims(data); % get the model dimensions

        muX = varargin{end-1};
        muY = varargin{end};
        
        % convert back into dimensions to apply boundary
        if isfield(data, 'targets')
            targs = data.targets;
        elseif isfield(data, 'items') && isfield(data, 'whichIsTestItem')
            whichItems = data.whichIsTestItem;
            for i = 1:size(data.items,2)
                targs(:,i) = data.items([whichItems(i)*2-1, whichItems(i)*2],i);
            end
        else
            error('need either targets or items & whichIsTestItem specified')
        end
        
        if isfield(data, 'errors')
            respLoc = data.errors + targs - [muX;muY];
            respLoc = bound(respLoc, dimensions);
            data.errors = respLoc - targs;
        end
        if isfield(data, 'changeSize')
            respLoc = data.changeSize + targs - [muX;muY];
            respLoc = bound(respLoc, dimensions);
            data.changeSize = respLoc - targs;
        end
        [p, pSep] = model.oldPdf(data, varargin{1:end-2});
    end
end

function t = bound(t,dimensions)
if any(t > dimensions | t < 0,'all')
    for i = 1:2
        t( i, t(i,:)>dimensions(i) ) = dimensions(i);
        t( i, t(i,:)<0 ) = 0;
    end
end

end
