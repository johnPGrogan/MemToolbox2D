% ENSEMBLEINTEGRATIONMODEL2D integration with distractors shifts reports
% Based on the model of Brady & Alvarez (2011), this model shifts
% observers' representations towards the mean values of the distractors,
% with more shift occuring if the distractors are closer together/closer to
% the target.
%
% In addition to data.errors, the data struct should include:
%   data.distractors, Row 1: distance of distractor 1 from target
%   ...
%   data.distractors, Row N: distance of distractor N from target
%
% Note that these values are the *distance from the correct answer*, not
% the actual color values of the distractors. They should thus range from
% -maxScreenDistance to maxScreenDistance.
%
% data.distractors may contain NaNs. For example, if you have data with
% different set sizes, data.distractors should contain as many rows as you
% need for the largest set size, and for displays with smaller set sizes
% the last several rows can be filled with NaNs.
%
% This is a simplified version of the Brady & Alvarez (2011) model, in that
% it does not take into account noise in the sampling of the distractors,
% and allows only a single level of integration (with all distractors).
%
% This 2D version is adapted from code in the MemToolbox package(Suchow, J. W.,
% Brady, T. F., Fougnie, D., & Alvarez, G. A. (2013). Modeling visual working 
% memory with the MemToolbox. Journal of Vision, 13(10):9, 1–8. 
% doi:10.1167/13.10.9. MemToolbox.org). 
% It was adapted by John P. Grogan, 2019.
%

function model = EnsembleIntegrationModel2D(boundary, maxIter)
  model.name = 'Ensemble integration model 2D';
	model.paramNames = {'g', 'sd', 'samples'};
	model.lowerbound = [0 0 0]; % Lower bounds for the parameters
	model.upperbound = [1 Inf Inf]; % Upper bounds for the parameters
	model.movestd = [0.02, 0.1, 0.05];
	model.pdf = @IntegrationModel2DPDF;
	model.start = [0.2, 10, 2;  % g, B, sd
    0.4, 15, 1;  % g, B, sd
    0.1, 20, 5]; % g, B, sd
	model.generator = @IntegrationModel2DGenerator;

    if exist('boundary','var')
        %   optional argument to use a hard boundary for generative samples outside
        %   screen, otherwise they are resampled by default
        model.boundary = boundary;
    end
    if exist('maxIter','var')
        %         maximum number of iterations if sampling, default = 100
        model.maxIter = maxIter;
    end
    
  % To specify a prior probability distribution, change and uncomment
  % the following line, where p is a vector of parameter values, arranged
  % in the same order that they appear in model.paramNames:
  % model.prior = @(p) (1);


function [p, pSep] = IntegrationModel2DPDF(data, g, sd, samples)
	if(~isfield(data, 'distractors'))
	  error('The integration model requires that you specify the distractors.')
    end
    dimensions = GetModelDims(data); % get the model dimensions


	data.distractors(end+1:end+2, :) = 0; % (target is always at zero)
	ensembleMean = [nanmean(data.distractors(1:2:end,:));nanmean(data.distractors(2:2:end,:))];
	ensembleStd = [nanstd(data.distractors(1:2:end,:));nanstd(data.distractors(2:2:end,:))];
	w = (1./(ensembleStd.^2)) ./ ((1./(ensembleStd.^2)) + (samples./(sd.^2))); % ratio of stimuli spread to samples over precision
	shiftedMean = w.*ensembleMean + (1-w)*0;  % (target is always at zero)

	pA = (1-g).*mvnpdf(data.errors',shiftedMean',[sd^2,sd^2]); % prob of target
	pG = (g).* (ones(length(data.errors),1) ./(prod(dimensions))); % prob of guess

	pSep = [pA, pG];

	p = sum(pSep, 2);
end

function r = IntegrationModel2DGenerator(parameters, dims, displayInfo)
% for simulating responses

  n = prod(dims); % figure out how many numbers to cook
  
  dimensions = GetModelDims(displayInfo); % get the model dimensions

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
  
  r = [rand(n,2) .* dimensions'] - targs'; % fill array with blind guesses
  
  guesses = logical(rand(n,1) < parameters{1}); % figure out which ones will be guesses

  displayInfo.distractors(end+1:end+2, :) = 0; % (target is always at zero)
  ensembleMean = [nanmean(displayInfo.distractors(1:2:end,:));nanmean(displayInfo.distractors(2:2:end,:))];
  ensembleStd = [nanstd(displayInfo.distractors(1:2:end,:));nanstd(displayInfo.distractors(2:2:end,:))];
  w = (1./(ensembleStd.^2)) ./ ((1./(ensembleStd.^2)) + (parameters{3}./(parameters{2}.^2))); % ratio of stimuli spread to samples over precision
  shiftedMean = w.*ensembleMean + (1-w)*0;  % (target is always at zero)

  
%   r(~guesses,:) = mvnrnd(shiftedMean(:,~guesses)',[parameters{2}^2,parameters{2}^2], sum(~guesses,1)); % target + noise
  
  % target errors - as response coords
  targResps = mvnrnd(targs' - shiftedMean', [parameters{2}^2,parameters{2}^2], n); % generate target response locations
  
  toReplace = any(targResps < 0 | targResps > dimensions',2); %check if out of bounds
  if any(toReplace,'all') %either resample (default) or apply boundary, depending on displayInfo fields
      targResps(toReplace,:) = ResampleOrBound(targResps(toReplace,:), targs(:,toReplace) - shiftedMean(:,toReplace), [parameters{2}^2,parameters{2}^2], model, dimensions);
  end
  r(~guesses,:) = targResps(~guesses,:)  - targs(:,~guesses)'; % targresp location - target = target error

  r = r';%reshape(r, dims); % reshape to requested dimensions
end

end