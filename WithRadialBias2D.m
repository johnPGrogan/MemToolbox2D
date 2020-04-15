% WITHRadialBIAS2D adds a radial bias towards any coordinate on each trial
% term to any 2D model. 
% data.biasCoords can be:
%   one coordinate [x;y] which is applied to all trials,
%   a vector of different coordiantes for each trial (i.e. [x(1:n); y(1:n)]
%   a number-NaN pair, so that the bias only applies to either X or Y
%   coordinates (e.g. [x; NaN] for a horizontal-only bias)
%   not given - the model will assume a bias towards the centre of the
%   screen
% 
% The bias parameter (mu) will be proportion of distance covered towards
% the biasCoords (positive values), or away from (negative values.
%
%  model = WithRadialBias2D(model, [priorForMu])
%
% The first parameter is the model to convert; the second (optional)
% parameter is a function for the prior on mu. It takes a single argument,
% mu.
%
% Thus, to take the StandardMixtureModel2D, which has a guess rate (g) and
% standard deviation (sd), and add a shift term (mu), just use:
%   model = WithRadialBias2D(StandardMixtureModel2D())
%
% This wrapper is compatible with both WithResponseSampling2D() and TwoAFC2D(). For
% example, the following works fine:
%   model = TwoAFC2D(WithRadialBias2D(StandardMixtureModel2D));
%
% This 2D version is adapted from code in the MemToolbox package(Suchow, J. W.,
% Brady, T. F., Fougnie, D., & Alvarez, G. A. (2013). Modeling visual working 
% memory with the MemToolbox. Journal of Vision, 13(10):9, 1–8. 
% doi:10.1167/13.10.9. MemToolbox.org). 
% It was adapted by John P. Grogan, 2019.
%

function model = WithRadialBias2D(model, priorForMu)
% If no prior is specified, default to an improper uniform prior
if nargin < 2
    priorForMu = @(p)(1);
end

% check that the model is a 2D model
if ~regexp('2D',model.name)
    warning('you are using WithBias2D() on a model which does not seem to be a 2D model')
end


% Take model and turn it into a model with a bias term
model.name = [model.name ' with radial bias 2D to coords'];
model.paramNames = [model.paramNames 'mu'];
model.lowerbound = [model.lowerbound -1];
model.upperbound = [model.upperbound 1];
model.movestd = [model.movestd 0.02];
model.start = [model.start rand(size(model.start,1),1) - .5];

% Adjust pdf and prior
model.oldPdf = model.pdf;
model.pdf = @NewPDF;

model.priorForMu = priorForMu;
if isfield(model, 'prior')
    model.oldPrior = model.prior;
    model.prior = @(p)(model.oldPrior(p(1:end-1)) .* model.priorForMu(p(end)));
end

% Adjust generator function
if isfield(model, 'generator')
    model.oldGenerator = model.generator;
    model.generator = @radialGenerator;%@(params,dims,displayInfo)(...
    %bound(model.oldGenerator(params(2:end), dims, displayInfo)+params{1}));
end

    function y = radialGenerator(params, dims, displayInfo)

        y = model.oldGenerator(params(1:end-1), dims, displayInfo); % get old generator
        
        y = applyRadialBias(y, displayInfo, params{end}); % apply the bias
    end

% Adjust model_plot
if isfield(model, 'modelPlot')
    model.oldModelPlot = model.modelPlot;
    model.modelPlot = @NewModelPlot;
end
    function figHand = NewModelPlot(data, params, varargin)
        if isstruct(params) && isfield(params, 'vals')
            mu = mean(params.vals(:,end));
            params.vals = params.vals(:, 1:end-1);
        else
            mu = params(end);
            params = params(1:end-1);
        end
        if isfield(data, 'errors')
            data.errors = applyRadialBias(data.errors, data, params); % apply the bias
        end
        if isfield(data, 'changeSize')
            data.changeSize = applyRadialBias(data.errors, data, params); % apply the bias;
        end
        figHand =  model.oldModelPlot(data, params, varargin);
    end

% Shift errors and/or changeSize
    function [p, pSep] = NewPDF(data, varargin)

        mu = varargin{end}; % get current mu
        mu = -mu/(1 - mu); % need to invert it to correct the radial bias
        if isfield(data, 'errors')
            data.errors = applyRadialBias(data.errors, data, mu); % apply the corrective bias;
        end
        if isfield(data, 'changeSize')
            data.changeSize = applyRadialBias(data.errors, data, mu); % apply the corrective bias
        end
        [p, pSep] = model.oldPdf(data, varargin{1:end-1});
    end
end

function t = bound(t,dimensions)
% t should be 2xnTrials
if any(t > dimensions | t < 0,'all')
    for i = 1:2
        t( i, t(i,:)>dimensions(i) ) = dimensions(i);
        t( i, t(i,:)<0 ) = 0;
    end
end
end

function err2 = applyRadialBias(err, displayInfo, mu)
% apply radial bias
dimensions = GetModelDims(displayInfo); % get the model dimensions


% need to know where the target was on each trial
if isfield(displayInfo, 'targets')
    targs = displayInfo.targets;
else
    if isfield(displayInfo,'items') && isfield(displayInfo,'whichIsTestItem') % if simulated data format
        items = displayInfo.items; % location of each item on each trial
        whichIs = displayInfo.whichIsTestItem; % which is chosen
        nItems = size(displayInfo.items,1)/2;
        for k=1:nItems % for each item
            targs(:,whichIs==k) = items([k*2-1,k*2],whichIs==k); % get location of target item
        end
    else
        error('cannot find target location, should be either data.targets or data.items and data.whichIsTestItem')
    end
end

% biasCoords
if isfield(displayInfo,'biasCoords')
    biasCoords = displayInfo.biasCoords; % get the coords to bias towards
else
    biasCoords = dimensions./2; % if none provided, fit bias to screen centre
end
if size(biasCoords,2) == 1 % if scalar, repmat, otherwise vector of biases for each trial
    biasCoords = repmat(biasCoords,1,size(err,2));
end

resp = targs + err;
diffs = biasCoords - resp; % get distance from response to coord bias

diffs(isnan(diffs)) = 0; % set any NaN coords to zero - wont' affect move

move = diffs .* mu; % get radial bias

resp2 = resp + move;
resp2 = bound(resp2, dimensions);

err2 = resp2 - targs;

% err2 = - targs + (resp + move); % apply bias
% err2 = bound(err2, dimensions); % reapply boundary

end