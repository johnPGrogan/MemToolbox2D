% TESTALLMODELSTWOAFC2D runs tests to be sure we are correctly sampling from and recovering the data
% for all the default 2D models when they are converted to 2AFC format
%
%    TestAllModelsTwoAFC2D(numTrials, numItemsPerTrial);
% 
% This 2D version is adapted from code in the MemToolbox package(Suchow, J. W.,
% Brady, T. F., Fougnie, D., & Alvarez, G. A. (2013). Modeling visual working 
% memory with the MemToolbox. Journal of Vision, 13(10):9, 1–8. 
% doi:10.1167/13.10.9. MemToolbox.org). 
% It was adapted by John P. Grogan, 2019.
%

function TestAllModelsTwoAFC2D(numTrials, numItemsPerTrial)

  % Default parameters
  if nargin < 1
    numTrials = 1000;
  end
  if nargin < 2
    numItemsPerTrial = 3;
  end

  % Which models to check:
  models = {...
    TwoAFC2D(StandardMixtureModel2D()), ...
    TwoAFC2D(WithBias2D(StandardMixtureModel2D)), ...
    TwoAFC2D(WithRadialBias2D(StandardMixtureModel2D)),...
    TwoAFC2D(SwapModel2D()), ...
  };

  % Try recovering parameters for each model
  for md=1:length(models)
    fprintf('\nModel: %s\n', models{md}.name);
    for s=1:size(models{md}.start,1)
      fprintf(' -- trying to recover params (%d of %d): ', s, size(models{md}.start,1));

      % Try sampling and fitting this model at each value of its start
      % parameters:
      paramsIn = models{md}.start(s,:);
      asCell = num2cell(models{md}.start(s,:));
      [paramsOut, lowerCI, upperCI] = ...
        TestSamplingAndFitting2D(models{md}, asCell, numTrials, ...
        numItemsPerTrial, 'Verbosity', 0);

      % Check that the credible intervals contain the correct parameter
      if all(paramsIn > lowerCI) && all(paramsIn < upperCI)
        fprintf('PASS\n');
      else
        which = find((paramsIn > lowerCI & paramsIn < upperCI) == 0, 1);
        fprintf('FAIL: %0.2f not in <%0.2f, %0.2f>\n', paramsIn(which), ...
          lowerCI(which), upperCI(which));
      end
    end
  end
end

