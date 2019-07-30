% SAMPLEFROMMODEL2D simulates data from a model with some parameters.
%
%  samp = SampleFromModel2D(model, params, dims, displayInfo)
%
%  Example usage:
%
%    model = StandardMixtureModel2D();
%    paramsIn = {0, 1};
%    simulatedData = SampleFromModel2D(model, paramsIn, [1,1000]);
%    paramsOut = MCMC(simulatedData, model);
%
%  If the model you pass in requires extra information for its pdf, for
%  example SwapModel2D (which requires .distractors), then you must pass in
%  display information as the 4th parameter (e.g., a data struct that has
%  .distractors), and the dimensions you request back (third parameter,
%  dims) must match the number of displays you provide.
%
% e.g.,
%    model = SwapModel2D();
%    displays = GenerateDisplays2D(100, 3); % 100 trials, 3 items/trial
%    paramsIn = {0, 0.1, 1};
%    simulatedData = SampleFromModel2D(model, paramsIn, [100,1], displays);
% 
% This 2D version is adapted from code in the MemToolbox package(Suchow, J. W.,
% Brady, T. F., Fougnie, D., & Alvarez, G. A. (2013). Modeling visual working 
% memory with the MemToolbox. Journal of Vision, 13(10):9, 1–8. 
% doi:10.1167/13.10.9. MemToolbox.org). 
% It was adapted by John P. Grogan, 2019.
%

function samp = SampleFromModel2D(model, params, dims, displayInfo)
  % Default to 1 sample
  if(nargin < 3)
    dims = [1 1];
  end

  % Make sure params is a cell array of parameters
  if ~iscell(params)
    params = num2cell(params);
  end

  % Check if the model needs extra information about the displays
  r = DoesModelRequireExtraInfo2D(model);
  if nargin < 4
    displayInfo = [];
    if r
      error(['You passed a model that requires extra information to make ' ...
        'a pdf; for example, maybe the set size (data.n) or the distractor ' ...
        'locations (data.distractors). Please pass the data struct describing ' ...
        'the displays you wish to sample data for as the fourth parameter.']);
    end
  end

  if r
    if isfield(displayInfo, 'errors')
      sz = size(displayInfo.errors,2);
    elseif isfield(displayInfo, 'distractors')
      sz = [1 size(displayInfo.distractors,2)];
    elseif isfield(displayInfo, 'n')
      sz = size(displayInfo.n);
    elseif isfield(displayInfo, 'afcCorrect')
      sz = size(displayInfo.afcCorrect);
    end
    if all(prod(dims) ~= prod(sz))
      error(['You passed a model that requires extra information to make ' ...
      'a pdf. When you pass such a model, dims needs to match the number of ' ...
      'displays you provide in the fourth parameter.']);
    end
  end

  % If the model has an efficient generator, use it. otherwise use rejection sampling
  if(isfield(model, 'generator'))
    samp = model.generator(params, dims, displayInfo);
    return
  end
  
    dimensions = GetModelDims(displayInfo); % get screen dims
    
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
    
    for j=1:2
        interpCoords(j,:) = linspace(0, dimensions(j), 100); % coords across screen
    end
    % combine all combinations of interpCoords - entire space
    interpVals = [reshape(repmat(interpCoords(1,:), 1, 100),[],1), reshape(repmat(interpCoords(2,:),100,1),[],1)]';

    data = displayInfo; % make a copy

    for i=1:length(interpVals) % for each point in space
      data.errors = repmat(interpVals(:,i),  dims) - data.targets; % get that distance form each trial's target
      y(i,:) = model.pdf(data, params{:}); % get pdf of that point for each trial
    end
    
    y2 = y ./ max(y); % normalise so max on each trial is 1
    y2 = permute(reshape(y2, 100, 100, []),[2,1,3]); % reshape
    
    u = rand(1,length(interpVals)); % random pdfvalues to take
    for i = 1:size(y2,3)
        % get coords that give a random pdfVal
        c = contourc(interpCoords(1,:), interpCoords(2,:), y2(:,:,i),[u(i) 1]);
        pdfVals = c(:,2:c(2,1)+1); % extract just those
        respLoc(:,i) = pdfVals(:, randi(size(pdfVals,2))); % pick one randomly
    end
    
    samples = respLoc - targs; % get distance of these responses from targets
    samp = reshape(samples, [], dims(2));
end


