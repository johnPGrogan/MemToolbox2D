% ALLGUESSINGMODEL2D returns a structure for a single-component model (a uniform
% distribution) for 2D data. This is like StandardMixtureModel2D, but without a
% remembered state.
%
% This 2D version is adapted from code in the MemToolbox package(Suchow, J. W.,
% Brady, T. F., Fougnie, D., & Alvarez, G. A. (2013). Modeling visual working 
% memory with the MemToolbox. Journal of Vision, 13(10):9, 1–8. 
% doi:10.1167/13.10.9. MemToolbox.org). 
% It was adapted by John P. Grogan, 2019.
%

function model = AllGuessingModel2D()
  model.name = 'All guessing model 2D';
	model.paramNames = {};
	model.lowerbound = []; % Lower bounds for the parameters
	model.upperbound = []; % Upper bounds for the parameters
	model.movestd = [];
%	model.pdf = @(data) 1*unifpdf(data.errors(:),-180,180);
	model.pdf = @AllGuessingModel2DPDF;
	model.start = []; % g
%   model.generator = @(parameters,dims,displayInfo) (unifrnd(-180,180,dims));
	model.generator = @AllGuessingModel2DGenerator;
end

function p = AllGuessingModel2DPDF(data)

  dimensions = GetModelDims(data); % get the model dimensions


  p = 1.* (ones(length(data.errors),1) ./(prod(dimensions))); % prob of guess

end

function r = AllGuessingModel2DGenerator(parameters, dims, displayInfo)
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

  r = (rand(n,2) .* dimensions') - targs'; % random guesses across screen
  
  r = r';%reshape(r, dims); % reshape to requested dimensions
end

