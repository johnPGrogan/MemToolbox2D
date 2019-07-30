% WITHLAPSES2D adds a lapse ("inattentional") term to any model
%
%  model = WithLapses2D(model, [priorForLapseRate])
%
% The first parameter is the model to convert; the second (optional)
% parameter is a function for the prior on the lapse rate. It takes a single
% argument, lapseRate.
%
% The lapse parameter adds in the possibility that observers might
% randomly guess on some proportion of trials ('lapses'). This is
% functionally identical to the guess parameters that already exist in
% several models (e.g.,StandardMixtureModel2D).
%
% Thus, to take the NoGuessingModel2D, which has a standard deviation (sd),
% and add an additional lapse term (lapseRate), just use:
%
%   model = WithLapses2D(NoGuessingModel2D())
%
% This model would then be identical to the StandardMixtureModel2D (since
% it will consist of both an sd and guess parameter).
%
% This wrapper is compatible with both TwoAFC2D(). For
% example, the following works fine:
%
%   model = TwoAFC2D(WithLapses2D(NoGuessingModel2D()));
%
% This 2D version is adapted from code in the MemToolbox package(Suchow, J. W.,
% Brady, T. F., Fougnie, D., & Alvarez, G. A. (2013). Modeling visual working 
% memory with the MemToolbox. Journal of Vision, 13(10):9, 1–8. 
% doi:10.1167/13.10.9. MemToolbox.org). 
% It was adapted by John P. Grogan, 2019.
%

function model = WithLapses2D(model, priorForLapseRate)
  % If no prior is specified, default to an improper uniform prior
  if nargin < 2
    priorForLapseRate = @(p)(1);
  end

  % Warn if the model already has a guess rate or lapse rate. Note that the built-in
  % models that fit capacities/guesses across set sizes (which is distinct from lapses)
  % call those parameters 'capacity' and not 'g', so this should only warn about
  % parameters that are truly just lapses:
  if any(strcmp(model.paramNames, 'g')) || any(strcmp(model.paramNames, 'lapse'))
    fprintf(['Warning: You are adding a lapse parameter to a model that already ' ...
      'has one ("g", guess rate or "lapse"). This is almost certainly not what you ' ...
      'want to do, since the two parameters will trade-off perfectly.\n\n']);
  end

  % Take model and turn it into a model with a bias term
  model.name = [model.name ' with lapses 2D'];
  model.paramNames = {'lapse', model.paramNames{:}};
  model.lowerbound = [0 model.lowerbound];
  model.upperbound = [1 model.upperbound];
  model.movestd = [0.02 model.movestd];
  model.start = [rand(size(model.start,1),1)  model.start];

  % Adjust pdf and prior
  model.oldPdf = model.pdf;
  model.pdf = @NewPDF;

  model.priorForLapseRate = priorForLapseRate;
  if isfield(model, 'prior')
    model.oldPrior = model.prior;
    model.prior = @(p)(model.oldPrior(p(2:end)) .* model.priorForLapseRate(p(1)));
  end

  % Shift errors and/or changeSize
  function [p, pSep] = NewPDF(data, lapseRate, varargin)
    dimensions = GetModelDims(data); % get the model dimensions


    p = (1-lapseRate) .* model.oldPdf(data, varargin{:});

    pG = (lapseRate) .* (ones(size(data.errors,2),1) ./(prod(dimensions))); % prob of guess
    
	pSep = [p, pG];
    p = sum(pSep,2);
  end
end
