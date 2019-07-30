% TWOAFC2D converts a 2D model so that it can be fit to 2AFC data
%
% Once a model has been converted, it no longer uses data.errors, but instead
% requires two new fields in the data struct:
%   data.afcCorrect - a sequence of zeros and ones, saying whether observers
%        got a 2AFC trial correct or incorrect
%   data.changeSize - [x;y] distance from correct to foil item in the 2FC
%        trials. Should correspond to .afcCorrect.
%
% TwoAFC2D can be wrapped around any 2D model. However, it is much slower 
% if wrapped around a model that requires additional fields of data, like 
% data.n or data.distractors (because then it must calculate a cdf for each
% datapoint separately).
%
% This 2D version is adapted from code in the MemToolbox package(Suchow, J. W.,
% Brady, T. F., Fougnie, D., & Alvarez, G. A. (2013). Modeling visual working 
% memory with the MemToolbox. Journal of Vision, 13(10):9, 1–8. 
% doi:10.1167/13.10.9. MemToolbox.org). 
% It was adapted by John P. Grogan, 2019.
%

function model = TwoAFC2D(model, samplesToApproxCDF)
  % How many samples of pdf to take to estimate the cdf with
  if nargin < 2
    samplesToApproxCDF = 100; % this gets squared to make 2D space
  end

  % Check if we need extra information to call the pdf
  model.requiresSeparateCDFs = DoesModelRequireExtraInfo2D(model);
  model.isTwoAFC = true;

  % Take model and turn it into a 2AFC-model
  model.name = ['2AFC2D ' model.name];
  model.oldPdf = model.pdf;
  model.pdf = @NewPDF;
  model.generator = @NewGenerator;

  % Make a generator of afcCorrect given changeSize
  function samples = NewGenerator(params, dims, displayInfo)
    displayInfo.afcCorrect = ones(1,size(displayInfo.changeSize,2));
    [p, thetas] = model.pdf(displayInfo, params{:});
    samples = binornd(1,thetas)';
  end

  % Convert pdf into a 2AFC pdf
  function [p,thetas] = NewPDF(data, varargin)
    
    dimensions = GetModelDims(data); % get the model dimensions

    
    nTr = size(data.changeSize,2);
    data.probe = data.targets + data.changeSize; % coords of probe
    
    % [x;y] coords to interpolate at, across entire screen error space
    for j = 1:2
        interpCoords(j,:) = linspace(0, dimensions(j), samplesToApproxCDF);
    end
    
    % combine all combinations of interpCoords - entire space
    interpVals = [reshape(repmat(interpCoords(1,:), 1,samplesToApproxCDF),[],1), reshape(repmat(interpCoords(2,:),samplesToApproxCDF,1),[],1)]';
    
    %find all points that are closer to zero than to data.changeSize
    leftPt = data.changeSize./2 + data.targets; % halfway point
    grad = data.changeSize(2,:) ./ data.changeSize(1,:); % gradient of line to data.changesize
    intercept = leftPt(2,:) - leftPt(1,:).* (-1./grad); % where it crosses y axis
    
    for i = 1:nTr
        interpIndsPerTrial(i,:) = vecnorm(interpVals - data.targets(:,i)) < vecnorm(interpVals - data.probe(:,i));
    end
    
    
    if ~model.requiresSeparateCDFs
      % Fast way
      % Get cdf for model
      y = NaN(samplesToApproxCDF^2,nTr);
      for i = 1:nTr
          data.errors = interpVals - data.targets(:,i);
          y(:,i) = model.oldPdf(data, varargin{:});
      end
           
      % normalise pdf to make units the same and sum to 1
      y = y ./ sum(y);
      
      y(~interpIndsPerTrial') = NaN; % set values closer to changeSize to NaN
            
      thetas = nansum(y)'; % sum for all the valid coords

    else
      % Slow way
      % Get separate pdf matrices for each data point
      y = NaN(samplesToApproxCDF^2,nTr);
      for i=1:samplesToApproxCDF^2
        data.errors = repmat(interpVals(:,i), 1, nTr) - data.targets;
        y(i,:) = model.oldPdf(data, varargin{:});
      end
      
      % normalise pdf
      y = y ./ sum(y);
      
      y(~interpIndsPerTrial') = NaN; % remove coords closer to changeSize than zero
      
      thetas = nansum(y)'; % sum pdf of valid coords

    end

    thetas(thetas>1) = 1;
    thetas(thetas<0) = 0;
    p = binopdf(data.afcCorrect(:), 1, thetas(:));
  end
end
