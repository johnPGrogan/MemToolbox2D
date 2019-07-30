% SLOTMODEL2D returns a structure for a two-component mixture model
% capacity K and precision sd. Capacity is the maximum number of independent
% representations. If the set size is greater than capacity some guesses will
% occur. For example, if participants can store 3 items but have to remember 6,
% participants will guess 50% of the time. Precision is the uncertainty of
% stored representations, and is assumed to be constant across set size.
%
% In addition to data.errors, requires data.n (the set size for each trial)
%
% A prior probability distribution can be specified in model.prior. Example
% priors are available in MemModels/Priors.
%
% This 2D version is adapted from code in the MemToolbox package(Suchow, J. W.,
% Brady, T. F., Fougnie, D., & Alvarez, G. A. (2013). Modeling visual working 
% memory with the MemToolbox. Journal of Vision, 13(10):9, 1–8. 
% doi:10.1167/13.10.9. MemToolbox.org). 
% It was adapted by John P. Grogan, 2019.
%

function model = SlotModel2D()
  model.name = 'Slot model 2D';
	model.paramNames = {'capacity', 'sd'};
	model.lowerbound = [0 0];     % Lower bounds for the parameters
	model.upperbound = [Inf Inf]; % Upper bounds for the parameters
	model.movestd = [0.25, 0.1];
	model.pdf = @slotpdf2D;
	model.start = [1, 20;   % capacity, sd
                 4, 40;  % capacity, sd
                 6, 60]; % capacity, sd

  % To specify a prior probability distribution, change and uncomment
  % the following line, where p is a vector of parameter values, arranged
  % in the same order that they appear in model.paramNames:
  % model.prior = @(p) (1);
end

function [p, pSep] = slotpdf2D(data,capacity,sd)
  
  if ~isfield(data, 'n')
      warning('no data.n field is supplied. Estimating set size from data.distractors')
      data.n = 1 + sum(~isnan(data.distractors),1)./2;
  end
  
  dimensions = GetModelDims(data); % get the model dimensions

  
  g = (1 - max(0,min(1,capacity./data.n(:))));

  pA = (1-g).*mvnpdf(data.errors',[0,0],[sd^2,sd^2]); % prob of target
  pG = (g).* (ones(size(data.errors,2),1) ./(prod(dimensions))); % prob of guess

  pSep = [pA, pG];
  
  p = sum(pSep,2);
end
