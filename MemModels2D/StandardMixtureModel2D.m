% STANDARDMIXTUREMODEL2D returns a structure for a two-component mixture model
% with guess rate g and standard deviation sd.
% for fitting 2D data (data.errors = [x;y])
%
% This version works on 2D data (e.g. x and y co-ordinates), and needs
% data.errors to be 2xN size, where the first row is X coords, and 2nd is Y
% coords.
%
%
% This version works on 2D data (e.g. x and y co-ordinates), and needs
% data.errors to be 2xN size, where the first row is X coords, and 2nd is Y
% coords.% This 2D version is adapted from code in the MemToolbox package(Suchow, J. W.,
% Brady, T. F., Fougnie, D., & Alvarez, G. A. (2013). Modeling visual working 
% memory with the MemToolbox. Journal of Vision, 13(10):9, 1–8. 
% doi:10.1167/13.10.9. MemToolbox.org). 
% It was adapted by John P. Grogan, 2019.
%

function model = StandardMixtureModel2D(boundary, maxIter)
    model.name = 'Standard mixture model 2D';
    model.paramNames = {'g', 'sd'};
    model.lowerbound = [0 0]; % Lower bounds for the parameters
    model.upperbound = [1 Inf]; % Upper bounds for the parameters
    model.movestd = [0.02, 0.1];
    model.pdf = @StandardMixtureModel2DPDF;

    model.start = [0.2, 10;  % g, sd
                 0.4, 15;  % g, sd
                 0.1, 20]; % g, sd
             
    model.generator = @StandardMixtureModel2DGenerator;
  
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



function [p, pSep] = StandardMixtureModel2DPDF(data, g, sd)
  % Parameter bounds check
  if g > 1
    p = zeros(size(data.errors,2),1);
    pSep = repmat(p,1,2);
    return;
  end

  dimensions = GetModelDims(data); % get the model dimensions
  
  pA = (1-g).*mvnpdf(data.errors',[0,0],[sd^2,sd^2]); % prob of target
  pG = (g).* (ones(length(data.errors),1) ./(prod(dimensions))); % prob of guess

  pSep = [pA,  pG];
  p = sum(pSep,2);
end

function y = StandardMixtureModel2DGenerator(parameters, dims, displayInfo)
% for simulating responses

  n = prod(dims); % figure out how many numbers to cook
  
  dimensions = GetModelDims(displayInfo); % get the model dimensions
  
  %get the targ location on each trial
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
  
  % Assign types to trials
  r = rand(n,1);
  whichType = zeros(n,1); % default = remembered target
  
  whichType(r<parameters{1}) = -1; % to make guess
  
  y = zeros(n,2); % Fill in with errors
  
  guessResps = rand(n, 2) .* dimensions'; % random guess locations
  y(whichType==-1,:) = guessResps(whichType==-1,:) - targs(:,whichType==-1)'; % random locations minus targets

  % target errors - as response coords
  targResps = mvnrnd(targs', [parameters{2}^2,parameters{2}^2], n); % generate target response locations

  toReplace = any(targResps < 0 | targResps > dimensions',2); %check if out of bounds
  if any(toReplace,'all') %either resample (default) or apply boundary, depending on displayInfo fields
      targResps(toReplace,:) = ResampleOrBound(targResps(toReplace,:), targs(:,toReplace), [parameters{2}^2,parameters{2}^2], model, dimensions);
  end
  y(whichType==0,:) = targResps(whichType==0,:)  - targs(:,whichType==0)'; % targresp location - target = target error
 
  y = y';%reshape(r, dims); % reshape to requested dimensions
end

end