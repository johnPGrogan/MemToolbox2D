% VARIABLEPRECISIONMODEL2D_GAMMASD returns a structure for a variable precision mixture model
% in which the standard deviations of observers' reports are assumed to be
% distributed as a gamma distribution.
%
% I've parameterized the gamma with a mode and SD, rather than the more
% traditional shape and rate, because this results in much better behaved
% posterior distributions and more interpretable parameters.
% 
% This 2D version is adapted from code in the MemToolbox package(Suchow, J. W.,
% Brady, T. F., Fougnie, D., & Alvarez, G. A. (2013). Modeling visual working 
% memory with the MemToolbox. Journal of Vision, 13(10):9, 1–8. 
% doi:10.1167/13.10.9. MemToolbox.org). 
% It was adapted by John P. Grogan, 2019.
%


function model = VariablePrecisionModel2D_GammaSD(boundary, maxIter)
  model.name = 'Variable precision model 2D (gamma over sd)';
	model.paramNames = {'g',  'modeSTD', 'sdSTD'};
	model.lowerbound = [0 0.01 0.01]; % Lower bounds for the parameters
	model.upperbound = [1 100 100]; % Upper bounds for the parameters
	model.movestd = [0.02, 0.15, 0.25];
	model.pdf = @vp_pdf;
    model.modelPlot = @model_plot;
	model.start = [0.0, 15, 5;
                 0.2, 20, 10;
                 0.1, 10, 2;
                 0.2, 30, 3];
    model.generator = @vp_gen;
  
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

  % For speed, calculate these all out here
  stdsSumOver = linspace(0.5, 100, 500);
  kValues = deg2k(stdsSumOver)';
  baseK = log(besseli(0, kValues, 1)) + kValues;
  lastX = [];

  function y = vp_pdf(data,g,modeSTD,sdSTD)
    dims = GetModelDims(data); % get the model dimensions
    
    % Probability of each of these
    scale = (2*sdSTD^2) / (modeSTD+sqrt(modeSTD^2+4*sdSTD^2));
    shape = 1 + modeSTD*(1/scale);
    probEachSD = gampdf(stdsSumOver, shape, scale);
    probEachSD = probEachSD./sum(probEachSD);

    if length(data.errors)~=length(lastX) || any(data.errors~=lastX,'all')
      % Calculate pdf for each STD; only do if the data is different than
      % last time
      model.v = [];
      x = repmat(data.errors', [1 1 length(stdsSumOver)]);
      k = repmat(stdsSumOver, [length(data.errors) 1]);
      for i=1:length(stdsSumOver) % targ pdf
          model.v(:,i) = mvnpdf(x(:,:,i), [0 0], [stdsSumOver(i).^2 stdsSumOver(i).^2]);
      end
      lastX = data.errors;
    end

    % Make final model prediction and sum
    probDataUnderThisNormal = (1-g).*model.v + (g).* (ones(size(model.v)) ./(prod(dims)));
    probEachSDBig = repmat(probEachSD, [size(probDataUnderThisNormal,1), 1]);
    y = sum(probDataUnderThisNormal.*probEachSDBig,2);
  end

  % Use our custom modelPlot to make a higher-order distribution plot
  function figHand = model_plot(data, params, varargin)
    figHand = figure();
    if isstruct(params) && isfield(params, 'vals')
      maxParams = MCMCSummarize(params, 'maxPosterior');
      params = params.vals(randsample(size(params.vals,1), 100),:);
    else
      maxParams = params;
    end
    set(gcf, 'Color', [1 1 1]);
    x = stdsSumOver;
    for i=1:size(params,1)
      modeSTD = params(i,2);
      sdSTD = params(i,3);
      scale = (2*sdSTD^2) / (modeSTD+sqrt(modeSTD^2+4*sdSTD^2));
      shape = 1 + modeSTD*(1/scale);
      y = gampdf(stdsSumOver, shape, scale);
      plot(x, y, 'Color', [0.54, 0.61, 0.06]); hold on;
    end
    modeSTD = maxParams(2);
    sdSTD = maxParams(3);
    scale = (2*sdSTD^2) / (modeSTD+sqrt(modeSTD^2+4*sdSTD^2));
    shape = 1 + modeSTD*(1/scale);
    y = gampdf(stdsSumOver, shape, scale);
    plot(x, y, 'k', 'LineWidth', 3);
    title('Higher-order distribution', 'FontSize', 14);
    xlabel('Standard dev. (degrees)', 'FontSize', 14);
    ylabel('Probability', 'FontSize', 14);
  end
    
  function [y, whichType] = vp_gen(params,dims,displayInfo)
    % Swap model random number generator

      if params{1} > 1
          error('params > 1')
      end
      n = prod(dims);
      dimensions = GetModelDims(displayInfo); % get the model dimensions
      
      numDistractorsPerTrial = [sum(~isnan(displayInfo.distractors),1)./2]';

      % distrib of stds
      scale = (2*params{3}^2) / (params{2} +sqrt(params{2}^2 + 4*params{3}^2));
      shape = 1 + params{2}*(1/scale);
      eachSD = gamrnd(shape, scale, n,1); % random stds form dist

      % Assign types to trials
      r = rand(n,1);
      whichType = zeros(n,1); % default = remembered target
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
      
      % target errors - as response coords
      sdVec = permute([eachSD.^2,eachSD.^2],[3,2,1]);
      targResps = mvnrnd(targs',  sdVec, n); % generate target response locations

      toReplace = any(targResps < 0 | targResps > dimensions',2); %check if out of bounds
      if any(toReplace,'all') %either resample (default) or apply boundary, depending on displayInfo fields
        targResps(toReplace,:) = ResampleOrBound(targResps(toReplace,:), targs(:,toReplace), sdVec(:,:,toReplace), model, dimensions);
      end
      y(whichType==0,:) = targResps(whichType==0,:)  - targs(:,whichType==0)'; % targresp location - target = target error
 
      y(whichType==-1,:) = (rand(sum(whichType==-1), 2) .* dimensions') ...
                            - targs(:,whichType==-1)'; % random locations minus targets

      % Reshape for output
      y = y';
    end
end

