% SLOTPLUSAVERAGINGMODEL returns a structure for the slots+averaging model of
% Zhang & Luck (2008). The model has two parameters: capacity, which is
% the number of available slots, and SD, which is the standard deviation
% of memory for an item that gets one slot.
%
% In addition to data.errors, requires data.n (the set size for each
% trial). The model is not particularly well-formed unless you have tested
% multiple set sizes; with only a single set size you may be better off
% with a model that does not make predictions across set size, like
% StandardMixtureModel2D().
%
% Uses the capacity and SD to fit data across multiple sizes.
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

function model = SlotsPlusAveragingModel2D()
  model.name = 'Slots+averaging model 2D';
	model.paramNames = {'capacity', 'sd'};
	model.lowerbound = [1 0];     % Lower bounds for the parameters
	model.upperbound = [Inf Inf]; % Upper bounds for the parameters
	model.movestd = [0.1, 0.1];
	model.pdf = @slotpdf2D;
	model.start = [2, 5;    % capacity, sd
                 3, 10;
                 4, 100];

  % To specify a prior probability distribution, change and uncomment
  % the following line, where p is a vector of parameter values, arranged
  % in the same order that they appear in model.paramNames:
  % model.prior = @(p) (1);
end

function [p, pSep] = slotpdf2D(data,capacity,sd)

  dimensions = GetModelDims(data); % get the model dimensions


  % First compute the number of items that get at least one slot
  numRepresented = min(capacity, data.n(:));

  % ... which we can use to compute the guess rate
  g = 1 - numRepresented ./ data.n(:);

  % Then pass around the slots evenly and compute the sd
  slotsPerItemEvenly = floor(capacity ./ data.n(:));
  worseSD = sd ./ max(sqrt(slotsPerItemEvenly),1); % to avoid Infs

  % Count the items that get an extra slot and the resulting sd
  numItemsWithExtraSlot = mod(capacity, data.n(:));
  pExtraSlot = numItemsWithExtraSlot ./ numRepresented;
  betterSD = sd ./ sqrt(slotsPerItemEvenly+1);

  % Finally, compute probability
  pA1 = (1-g(:)).*(pExtraSlot(:)) .*mvnpdf(data.errors',[0,0],permute([betterSD.^2, betterSD.^2], [3,2,1])); % prob of target
  pA2 = (1-g(:)).*(1-pExtraSlot(:)) .*mvnpdf(data.errors',[0,0],permute([worseSD.^2, worseSD.^2], [3,2,1])); % prob of target
  pG = (g).* (ones(size(data.errors,2),1) ./(prod(dimensions))); % prob of guess

  pSep = [pA1, pA2, pG];
  p = sum(pSep,2);
end
