% SWAPMODELSepMisbind2D 
% allows separate misbinding parameters for two different sets of
% non-target items, e.g. non-probed targets and distractors, or previous
% response/target and distractors.
% Based on Bays, Catalao, & Husain (2009) model.
% 
% These two sets of distance co-ordinates should be entered as
%   data.distractors1, Row 1: x coord distance of distractor1 1 from target
%                      Row 2: y coord distance of distractor1 1 from target
%   ...
%   data.distractors1, Row N-1: x coord distance of distractor1 N from target
%                      Row N  : y coord distance of distractor1 N from target
% 
% and the same for data.distractors2
% If you wish to fit misbinding to distractors and previous responses
% separately, you can set data.distractors1 = data.distractors, and
% data.distractors2 = data.previousResponse.
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
%
% This 2D version is adapted from code in the MemToolbox package(Suchow, J. W.,
% Brady, T. F., Fougnie, D., & Alvarez, G. A. (2013). Modeling visual working 
% memory with the MemToolbox. Journal of Vision, 13(10):9, 1–8. 
% doi:10.1167/13.10.9. MemToolbox.org). 
% It was adapted by John P. Grogan, 2019.
%

function model = SwapModelSepMisbind2D(boundary, maxIter)
    model.name = 'Swap model 2D with separate misbinding to 2 different sets of non-targets';
	model.paramNames = {'g', 'B1','B2' 'sd'};
	model.lowerbound = [0 0 0 0]; % Lower bounds for the parameters
	model.upperbound = [1 1 1 Inf]; % Upper bounds for the parameters
	model.movestd = [0.02, 0.02, 0.02, 0.5];
    model.pdf = @SwapModel2DPDF;
    model.modelPlot = @model_plot;
    model.generator = @SwapModel2DGenerator;
    model.start = [0.2, 0.1, 0.1, 20;  % g, BT, BD, sd
                   0.4, 0.1, 0.2, 40;  % g, BT, BD, sd
                   0.1, 0.5, 0.4, 60]; % g, BT, BD, sd

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
  end


function [p, pSep] = SwapModel2DPDF(data, g, B1, B2, sd)
  % Parameter bounds check
  if g+B1+B2 > 1
    p = zeros(size(data.errors,2),1);
    pSep = repmat(p, 1,4);
    return;
  end

  if(~isfield(data, 'distractors1')) || ~isfield(data,'distractors2')
    error('The swap model requires that you specify distractors1 and distractors2.')
  end
  
  dimensions = GetModelDims(data); % get the model dimensions


  pA = (1-g-B1-B2).*mvnpdf(data.errors',[0,0],[sd^2,sd^2]); % prob of target
  pG = (g).* (ones(size(data.errors,2),1) ./(prod(dimensions))); % prob of guess

  
  % misbinding to distractors1
  numDist1PerTrial = [sum(~isnan(data.distractors1),1)./2]';
  pBD1 = NaN(size(data.errors,2), max(numDist1PerTrial)); % preallocate in case there are no distractors
  for i=1:max(numDist1PerTrial)
    pdfOut = mvnpdf(data.errors', data.distractors1((i*2-1):i*2,:)', diag([sd^2,sd^2],0)); % prob of misbind to this distractor
    pdfOut(isnan(pdfOut)) = 0;
    pBD1(:,i) = (B1./numDist1PerTrial).*pdfOut;%need to do this for each dist
  end
  
  % misbinding to distractors2
  numDist2PerTrial = [sum(~isnan(data.distractors2),1)./2]';
  pBD2 = NaN(size(data.errors,2), max(numDist2PerTrial)); % preallocate in case there are no distractors
  for i=1:max(numDist2PerTrial)
    pdfOut = mvnpdf(data.errors', data.distractors2((i*2-1):i*2,:)', diag([sd^2,sd^2],0)); % prob of misbind to this distractor
    pdfOut(isnan(pdfOut)) = 0;
    pBD2(:,i) = (B2./numDist2PerTrial).*pdfOut;%need to do this for each dist
  end
  
  pSep = [pA, pG, pBD1, pBD2];
  
  p = sum(pSep,2); % sum across diff response types to get p per trial
  
end


function [y, whichType] = SwapModel2DGenerator(params,dims,displayInfo)
% Swap model random number generator

if params{1} + params{2} + params{3} > 1
    error('params > 1')
end

if ~isfield(displayInfo,'distractors1') || ~isfield(displayInfo, 'distractors2')
    error('distractors1 and distractors2 not supplied')
end

  g = params{1};
  BD1 = params{2};
  BD2 = params{3};
  sd = params{4};
  a = 1 - g - BD1 - BD2;
  n = prod(dims);

  dimensions = GetModelDims(displayInfo); % get the model dimensions
  
  numDist1PerTrial = [sum(~isnan(displayInfo.distractors1),1)./2]';
  numDist2PerTrial = [sum(~isnan(displayInfo.distractors2),1)./2]';
  
  % Assign types to trials
  r = rand(n,1);
  whichType = zeros(n,1); % default = remembered target
  
  whichType(r<(g+BD1+BD2)) = ceil(rand(sum(r<(g+BD1+BD2)), 1) ...
    .* numDist1PerTrial(r<(g+BD1+BD2))); % to swap to random distractor1
  whichType(r<(g+BD2)) = ceil(rand(sum(r<(g+BD2)), 1) ...
    .* numDist1PerTrial(r<(g+BD2)))+.5; % to swap to random distractor2

  whichType(r<g) = -1; % to make guess

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
  targResps = mvnrnd(targs', [sd^2, sd^2], n); % generate target response locations
  
  toReplace = any(targResps < 0 | targResps > dimensions',2); %check if out of bounds
  if any(toReplace,'all') %either resample (default) or apply boundary, depending on displayInfo fields
      targResps(toReplace,:) = ResampleOrBound(targResps(toReplace,:), targs(:,toReplace), [sd^2, sd^2], model, dimensions);
  end
  y(whichType==0,:) = targResps(whichType==0,:)  - targs(:,whichType==0)'; % targresp location - target = target error
  
 
  %distractor1 errors
    for d=1:(size(displayInfo.distractors1,1)/2)
        if any(whichType==d)
            distLoc = targs + displayInfo.distractors((d*2-1):d*2,:); % location of distractors
            distResps = mvnrnd(distLoc', [sd^2, sd^2], n);

            toReplace = any(distResps < 0 | distResps > dimensions',2); %check if out of bounds
            if any(toReplace,'all') %either resample (default) or apply boundary, depending on displayInfo fields
                distResps(toReplace,:) = ResampleOrBound(distResps(toReplace,:), distLoc(:,toReplace), [params{3}^2,params{3}^2], model, dimensions);
            end
            
            y(whichType==d,:) = distResps(whichType==d,:) - targs(:,whichType==d)';

        end
    end
    
    %distractor2 errors
    for d=1:(size(displayInfo.distractors2,1)/2)
        if any(whichType==d+.5)
            toReplace = true(n,1); % make n distresps initially
            distResps = NaN(n,2);
            while any(toReplace) % keep going until all are within screen
                distResps(toReplace,:) = mvnrnd(targs(:,toReplace)' + ...
                    displayInfo.distractors2((d*2-1):d*2,toReplace)', [sd^2,sd^2], sum(toReplace));
                toReplace = any(distResps < 0 | distResps > dimensions',2);
            end
            y(whichType==d+.5,:) = distResps(whichType==d+.5,:) - targs(:,whichType==d+.5)';
        end
    end
    
    


%   %distractor errors
%     for d=1:(size(displayInfo.distractors1,1)/2)
%         if any(whichType==d)
%             y(whichType==d,:) = mvnrnd(displayInfo.distractors1((d*2-1):d*2,(whichType==d))', ...
%               [sd^2,sd^2], sum(whichType==d)); % distractor location + noise
%         end
%     end
%     %distractor errors
%     for d=1:(size(displayInfo.distractors2,1)/2)
%         if any(whichType==(d+.5))
%             y(whichType==(d+.5),:) = mvnrnd(displayInfo.distractors2((d*2-1):d*2,(whichType==(d+.5)))', ...
%               [sd^2,sd^2], sum(whichType==(d+.5))); % distractor location + noise
%         end
%     end
    


  
  % Reshape for output
  y = y';
end

end