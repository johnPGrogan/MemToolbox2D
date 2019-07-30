% CONTINUOUSRESOURCEMODEL2D returns a structure for a continuous resource model,
% very much like the one from Bays & Husain (2008).
%
% Usage:
%   MemFit2D(data, ContinuousResourceModel2D);
%
% This 2D version is adapted from code in the MemToolbox package(Suchow, J. W.,
% Brady, T. F., Fougnie, D., & Alvarez, G. A. (2013). Modeling visual working 
% memory with the MemToolbox. Journal of Vision, 13(10):9, 1–8. 
% doi:10.1167/13.10.9. MemToolbox.org). 
% It was adapted by John P. Grogan, 2019.
%

function model = ContinuousResourceModel2D()
  model.name = 'Continuous resource model 2D';
  model.paramNames = {'lapse','k', 'bestSD'};
  model.lowerbound = [0 0 0]; % Lower bounds for the parameters
  model.upperbound = [1 10 50]; % Upper bounds for the parameters
  model.movestd = [0.02, 0.1, 1];
  model.pdf = @crpdf2D;
  model.start = [0.2, 0.1, 10;  % lapse, k, bestSD
                 0.4, 1, 2;
                 0.1, 10, 20];

  % To specify a prior probability distribution, change and uncomment
  % the following line, where p is a vector of parameter values, arranged
  % in the same order that they appear in model.paramNames:
  % model.prior = @(p) (1);
end

function [p, pSep] = crpdf2D(data,lapse,k,bestSD)

  dimensions = GetModelDims(data); % get the model dimensions


  propResources = 1./data.n(:);
  pMax = 1./(bestSD.^2);
  precision = (propResources .^ k) .* pMax;
  sd = sqrt(1./precision);

  pA = (1-lapse) .* mvnpdf(data.errors',[0,0],permute([sd.^2, sd.^2], [3,2,1])); % prob of target
  pG = (lapse) .* (ones(size(data.errors,2),1) ./(prod(dimensions))); % prob of guess
  
  pSep = [pA, pG];
  
  p = sum(pSep,2);
end
