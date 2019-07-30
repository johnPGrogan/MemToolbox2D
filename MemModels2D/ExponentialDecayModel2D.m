% EXPONENTIALDECAYMODEL2D returns a struct for a model where objects drop out
% of memory indepdendently at a constant rate over time, somewhat like the
% model proposed by Zhang & Luck (2009), though without an initial period of
% stability. (This is a pure death process over objects.)
%
% Parameters: tau, K, sd
%   Tau is the mean lifetime of an object (in ms), K is the capacity of
%   working memory, and sd is the precision with which items are remembered.
%
% In addition to data.errors, requires data.n (the set size for each trial)
% as well as data.time (the delay before test in each trial; in milliseconds).
%
% This 2D version is adapted from code in the MemToolbox package(Suchow, J. W.,
% Brady, T. F., Fougnie, D., & Alvarez, G. A. (2013). Modeling visual working 
% memory with the MemToolbox. Journal of Vision, 13(10):9, 1–8. 
% doi:10.1167/13.10.9. MemToolbox.org). 
% It was adapted by John P. Grogan, 2019.
%

function model = ExponentialDecayModel2D()
  model.name = 'Exponential decay model 2D';
	model.paramNames = {'tau', 'K', 'sd'};
	model.lowerbound = [0 0 0]; % Lower bounds for the parameters
	model.upperbound = [Inf Inf Inf]; % Upper bounds for the parameters
	model.movestd = [20, 1, 1];
	model.pdf = @sdpdf2D;

	model.start = [1000, 4, 12;  % tau, k, sd
                 2000, 2, 20;
                 10000, 6, 30];

  % To specify a prior probability distribution, change and uncomment
  % the following line, where p is a vector of parameter values, arranged
  % in the same order that they appear in model.paramNames:
  % model.prior = @(p) (1);
end

function [p, pSep] = sdpdf2D(data, tau, k, sd)

  if ~isfield(data,'time')
      error('data.time not supplied');
  end

  dimensions = GetModelDims(data); % get the model dimensions

  
  B = min(k, data.n); % maximum contribution of working memory

  % the probability of remembering is exponential in time
  p = B.*exp(data.time/-tau) ./ data.n;
  g = 1 - p; % the guess rate

  pA = (1-g(:)).*mvnpdf(data.errors',[0,0],[sd^2,sd^2]); % prob of target
  pG = (g(:)).* (ones(size(data.errors,2),1) ./(prod(dimensions))); % prob of guess

  pSep = [pA, pG];
  
  p = sum(p, 2);
end

