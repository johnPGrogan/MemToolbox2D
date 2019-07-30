% SWAPMODELXYVAR2D returns a structure for a four-component model
% with guesses and swaps, and separate X and Y variance.
% Based on Bays, Catalao, & Husain (2009) model.
% This is an extension of the StandardMixtureModel that allows for
% observers' misreporting incorrect items.
%
% This version works on 2D data (e.g. x and y co-ordinates), and needs
% data.errors to be 2xN size, where the first row is X coords, and 2nd is Y
% coords.
% 
% In addition to data.errors, the data struct should include:
%   data.distractors, Row 1: x coord distance of distractor 1 from target
%                     Row 2: y coord distance of distractor 1 from target
%   ...
%   data.distractors, Row N-1: x coord distance of distractor N from target
%                     Row N  : y coord distance of distractor N from target
%
% Note that these values are the *distance from the correct answer*, not
% the actual coord values of the distractors. They should thus range from
% -maxScreenDim to maxScreenDim.
%
% data.distractors may contain NaNs. For example, if you have data with
% different set sizes, data.distractors should contain twice as many rows
% ([x;y]) as you need for the largest set size, and for displays with 
% smaller set sizes the last several rows can be filled with NaNs.
%
% This model includes a custom .modelPlot function that is called by
% MemFit2D(). This function produces a plot of the distance of observers'
% reports from the distractors, rather than from the target, as in Bays,
% Catalao & Husain (2009), Figure 2B.
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

function model = SwapModelXYVar2D(boundary, maxIter)
    model.name = 'Swap model 2D with X and Y variances';
	model.paramNames = {'g', 'B', 'sdX', 'sdY'};
	model.lowerbound = [0 0 0 0]; % Lower bounds for the parameters
	model.upperbound = [1 1 Inf Inf]; % Upper bounds for the parameters
	model.movestd = [0.02, 0.02, 0.5, 0.5];
    model.pdf = @SwapModel2DPDF;
    model.modelPlot = @model_plot;
    model.generator = @SwapModel2DGenerator;
    model.start = [0.2, 0.1, 20, 10;  % g, B, sd
                   0.4, 0.1, 40, 50;  % g, B, sd
                   0.1, 0.5, 60, 60]; % g, B, sd

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

  % Use our custom modelPlot to make a plot of errors centered on
  % distractors (ala Bays, Catalao & Husain, 2009, Figure 2B)
  function figHand = model_plot(data, params, varargin)
   if isfield(data,'errors')
    d.errors = [];
    for i=1:length(data.errors)
      d.errors = [d.errors; circdist(data.errors(i), data.distractors(:,i))];
    end
    if isstruct(params) && isfield(params, 'vals')
      params = MCMCSummarize(params, 'maxPosterior');
    end
    m = StandardMixtureModel2D();
    f = [1-(params(2)/size(data.distractors,1)) params(3)];
    figHand = PlotModelFit2D(m, f, d, 'NewFigure', true, 'ShowNumbers', false);
    title('Error relative to distractor locations', 'FontSize', 14);
    topOfY = max(ylim);
    txt = sprintf('B: %.3f\nsd: %0.2f\n', params(2), params(3));
    text(180, topOfY-topOfY*0.05, txt, 'HorizontalAlignment', 'right');
   else
       warning('no data.errors field. not drawing model_plot');
   end
  end


function [p, pSep] = SwapModel2DPDF(data, g, B, sdX, sdY)
  % Parameter bounds check
  if g+B > 1
    p = zeros(size(data.errors,2),1);
    pSep = repmat(p,1,3);
    return;
  end

  if(~isfield(data, 'distractors'))
    error('The swap model requires that you specify the distractors.')
  end
  
  dimensions = GetModelDims(data); % get the model dimensions
  
  pA = (1-g-B).*mvnpdf(data.errors',[0,0],[sdX^2,sdY^2]); % prob of target
  pG = (g).* (ones(size(data.errors,2),1) ./(prod(dimensions))); % prob of guess

  
  % Allow for the possibility of NaN's in distractors, as in the case where
  % different trials have different set sizes
  numDistractorsPerTrial = [sum(~isnan(data.distractors),1)./2]';
  
  
  pB = NaN(size(data.errors,2), max(numDistractorsPerTrial)); % preallocate in case there are no distractors
  for i=1:max(numDistractorsPerTrial)
    pdfOut = mvnpdf(data.errors', data.distractors((i*2-1):i*2,:)', diag([sdX^2,sdY^2],0)); % prob of misbind to this distractor
    pdfOut(isnan(pdfOut)) = 0;
    pB(:,i) = (B./numDistractorsPerTrial).*pdfOut;%need to do this for each dist
  end
  
  pSep = [pA, pG, pB];
  
  p = sum(pSep,2); % sum across diff response types to get p per trial
  
end


function [y, whichType] = SwapModel2DGenerator(params,dims,displayInfo)
% Swap model random number generator

if params{1} + params{2} > 1
    error('params > 1')
end

  n = prod(dims);

  dimensions = GetModelDims(displayInfo); % get the model dimensions
  
  numDistractorsPerTrial = [sum(~isnan(displayInfo.distractors),1)./2]';
  
  % Assign types to trials
  r = rand(n,1);
  whichType = zeros(n,1); % default = remembered target
  
  whichType(r<params{1}+params{2}) = ceil(rand(sum(r<params{1}+params{2}), 1) ...
    .* numDistractorsPerTrial(r<params{1}+params{2})); % to swap to random distractor
  whichType(r<params{1}) = -1; % to make guess

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
  
  y = zeros(n,2); % Fill in with errors
   
  guessResps = rand(n, 2) .* dimensions'; % generate random guess locations
  y(whichType==-1,:) = guessResps(whichType==-1,:) - targs(:,whichType==-1)'; % random locations minus targets
  
  % target errors - as response coords
  targResps = mvnrnd(targs', [params{3}^2,params{4}^2], n); % generate target response locations
  
  toReplace = any(targResps < 0 | targResps > dimensions',2); %check if out of bounds
  if any(toReplace,'all') %either resample (default) or apply boundary, depending on displayInfo fields
      targResps(toReplace,:) = ResampleOrBound(targResps(toReplace,:), targs(:,toReplace), [params{3}^2,params{4}^2], model, dimensions);
  end
  y(whichType==0,:) = targResps(whichType==0,:)  - targs(:,whichType==0)'; % targresp location - target = target error
  
  

  %distractor errors
    for d=1:(size(displayInfo.distractors,1)/2)
        if any(whichType==d)
            distLoc = targs + displayInfo.distractors((d*2-1):d*2,:); % location of distractors
            distResps = mvnrnd(distLoc', [params{3}^2,params{4}^2], n);

            toReplace = any(distResps < 0 | distResps > dimensions',2); %check if out of bounds
            if any(toReplace,'all') %either resample (default) or apply boundary, depending on displayInfo fields
                distResps(toReplace,:) = ResampleOrBound(distResps(toReplace,:), distLoc(:,toReplace), [params{3}^2,params{4}^2], model, dimensions);
            end
        
            y(whichType==d,:) = distResps(whichType==d,:) - targs(:,whichType==d)';
        end
    end


  
  % Reshape for output
  y = y';
end
end