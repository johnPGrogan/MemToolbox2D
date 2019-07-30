% SLOTSPLUSRESOURCESMODEL2D returns a structure for the slots+resources model of
% Zhang & Luck (2008), though it assumes even allocation of the resource to
% whichever objects are assigned a slot. The model has two parameters: capacity,
% which is the number of available slots, and bestSD, which is the standard
% deviation of memory for an item when all of the resource is thrown at it.
% model with capacity K and precision sd.
%
% In addition to data.errors, this requires data.n (the set size for each
% trial). The model is not particularly well-formed unless you have tested
% multiple set sizes; with only a single set size you may be better off
% with a model that does not make predictions across set size, like
% StandardMixtureModel2D().
%
% Uses the capacity and bestSD to fit data across multiple sizes.
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

function model = SlotsPlusResourcesModel2D()
  model.name = 'Slot plus resouces model 2D';
	model.paramNames = {'capacity', 'bestSD'};
	model.lowerbound = [0 0]; % Lower bounds for the parameters
	model.upperbound = [Inf Inf]; % Upper bounds for the parameters
	model.movestd = [0.1, 0.1];
	model.pdf = @slotpdf2D;
	model.start = [2, 20;  % g, sd
                   3, 60;  % g, sd
                   4, 5]; % g, sd

  % To specify a prior probability distribution, change and uncomment
  % the following line, where p is a vector of parameter values, arranged
  % in the same order that they appear in model.paramNames:
  % model.prior = @(p) (1);
end

function [p, pSep] = slotpdf2D(data,capacity,bestSD)

  dimensions = GetModelDims(data); % get the model dimensions


  numRepresented = min(capacity, data.n(:));
  g = 1 - numRepresented ./ data.n(:);
  sd = bestSD .* sqrt(numRepresented(:));

  pA = (1-g).*mvnpdf(data.errors',[0,0],permute([sd.^2, sd.^2], [3,2,1])); % prob of target
  pG = (g).* (ones(size(data.errors,2),1) ./(prod(dimensions))); % prob of guess

  pSep = [pA, pG];
  p = sum(pSep, 2);

end
