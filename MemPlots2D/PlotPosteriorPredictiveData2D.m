%PLOTPOSTERIORPREDICTIVEDATA2D Show data sampled from the model with the actual
% data overlayed, plus a plot of where the two differ. This can be thought
% of as the 'residual' of the data, given the model fit, and is helpful for
% diagnosing bad fits in the model.
%
%   figHand = PlotPosteriorPredictiveData2D(model, posteriorSamples,...
%                                             data, [optionalParameters])
%
%
% Optional parameters:
%  'NumberOfBins' - the number of bins to use in display the data. Default
%  55.
%
%  'NumSamplesToPlot' - how many posterior samples to show in the posterior
%  predictive plot. Default is 48.
%
%  'PdfColor' - the color to plot the model fit with.
%
%  'NewFigure' - whether to make a new figure or plot into the currently
%  active subplot. Default is false (e.g., plot into current plot).
% 
% This 2D version is adapted from code in the MemToolbox package(Suchow, J. W.,
% Brady, T. F., Fougnie, D., & Alvarez, G. A. (2013). Modeling visual working 
% memory with the MemToolbox. Journal of Vision, 13(10):9, 1–8. 
% doi:10.1167/13.10.9. MemToolbox.org). 
% It was adapted by John P. Grogan, 2019.
%

function figHand = PlotPosteriorPredictiveData2D(model, posteriorSamples, data, varargin)
  % Show data sampled from the model with the actual data overlayed, plus a
  % difference plot.
  args = struct('NumSamplesToPlot', 48, 'NumberOfBins', 55, ...
    'PdfColor', [0.54, 0.61, 0.06], 'NewFigure', true);
  args = parseargs(varargin, args);
  if args.NewFigure, figHand = figure(); else figHand = []; end

  % Choose which samples to use
  if(isempty(model.paramNames))
    which = 1:args.NumSamplesToPlot;
  else
    which = round(linspace(1, size(posteriorSamples.vals,1), args.NumSamplesToPlot));
  end

  dims = GetModelDims(data); % get the model dimensions


  
  % How to bin
  x = [linspace(-dims(1)/2, dims(1)/2, args.NumberOfBins)',linspace(-dims(2)/2, dims(2)/2, args.NumberOfBins)'];
  [normalizedData, nSamples] = getNormalizedBinnedData2D(data, x);

  for iDim=1:2
      % Plot samples
      subplot(2,2,iDim);
      set(gcf, 'Color', [1 1 1]);
      curFigure = gcf;
      hold on;
      sampTime = tic();
      curHandle = [];
      for i=1:length(which)
        % Generate random data from this distribution with these parameters
        if(isempty(model.paramNames))
          asCell = {};
        else
          asCell = num2cell(posteriorSamples.vals(which(i),:));
        end
        yrep = SampleFromModel2D(model, asCell, [1 nSamples], data);
        if i==1 && toc(sampTime)>(5.0/length(which)) % if it will take more than 5s...
          curHandle = awaitbar(i/length(which), ...
            'Sampling to get posterior predictive distribution...');
        elseif ~isempty(curHandle)
          if awaitbar(i/length(which))
            break;
          end
          set(0, 'CurrentFigure', curFigure);
        end

        % Bin data and model
        normalizedYRep = getNormalizedBinnedReplication(yrep, data, x);
        if any(isnan(normalizedYRep))
          hSim = plot(x(:,iDim), normalizedYRep(iDim,:), 'x-', 'Color', args.PdfColor);
          if verLessThan('matlab','8.4.0')
            set(hSim, 'LineSmoothing', 'on');
          end
        else
          hSim = patchline(x(:,iDim), normalizedYRep(iDim,:), 'LineStyle', '-', 'EdgeColor', ...
            args.PdfColor, 'EdgeAlpha', 0.15);
          if verLessThan('matlab','8.4.0')
            set(hSim, 'LineSmoothing', 'on');
          end
        end

        % Difference between this data and real data
        diffPlot((i*2-1):i*2,:) = normalizedData - normalizedYRep;
      end
      if ishandle(curHandle)
        close(curHandle);
      end

      % Plot data
      hSim = plot(-191:-190, [1 1], '-', 'Color', args.PdfColor);
      h=plot3(x(:,iDim),normalizedData(iDim,:),ones(size(x)),'ok-','LineWidth', 1.5, 'MarkerEdgeColor',[0 0 0], ...
           'MarkerFaceColor', [0 0 0], 'MarkerSize', 4);
      title('Simulated data from model');
      legend([h(1), hSim], {'Actual data', 'Simulated data'});
      legend boxoff;


      % Plot difference
      subplot(2,2,iDim+2);
      bounds = quantile(diffPlot(iDim:2:end,:), [.05 .50 .95])';
      if any(isnan(bounds))
        hB = errorbar(x(:,iDim), bounds(:,2), bounds(:,2)-bounds(:,1), bounds(:,3)-bounds(:,2), ...
          'x', 'Color', [.3 .3 .3], 'LineWidth', 2, 'MarkerSize', 10);
      else
        hB = boundedline(x(:,iDim), bounds(:,2), [bounds(:,2)-bounds(:,1) bounds(:,3)-bounds(:,2)], ...
          'cmap', [0.3 0.3 0.3]);
        set(hB, 'LineWidth', 2);
        if verLessThan('matlab','8.4.0')
          set(hB, 'LineSmoothing', 'on');
        end
      end
      line([-180 180], [0 0], 'LineStyle', '--', 'Color', [.5 .5 .5]);


      title('Difference between actual and simulated data');
      xlabel('(Note: deviations from zero indicate bad fit)');

  end
  
  % Allow the user to limit this figure to any subset of the data
  if ~isempty(figHand)
    CreateMenus(data, @redrawFig);
  end
  function redrawFig(whichField, whichValue)
    if strcmp(whichField, 'all')
      subplot(1,1,1);
      PlotPosteriorPredictiveData(model, posteriorSamples, data, ...
        'NewFigure', false, 'NumSamplesToPlot', args.NumSamplesToPlot, ...
        'NumberOfBins', args.NumberOfBins, 'PdfColor', args.PdfColor);
    elseif sum(ismember(data.(whichField),whichValue)) > 0
      [datasets,conditionOrder] = SplitDataByField(data, whichField);
      newData = datasets{ismember(conditionOrder,whichValue)};
      subplot(1,1,1);
      PlotPosteriorPredictiveData(model, posteriorSamples, newData, ...
        'NewFigure', false, 'NumSamplesToPlot', args.NumSamplesToPlot, ...
        'NumberOfBins', args.NumberOfBins, 'PdfColor', args.PdfColor);
    end
  end
end

function y = getNormalizedBinnedReplication(yrep, data, x)
for iDim=1:2
  if isfield(data, 'errors')
    y(iDim,:) = hist(yrep(iDim,:), x(:,iDim))';
    y(iDim,:) = y(iDim,:) ./ sum(y(iDim,:),'all');
  else
    for i=1:length(x)
      distM(:,i) = (data.changeSize(iDim,:) - x(i,iDim)).^2;
    end
    [tmp, whichBin] = min(distM,[],2);
    for i=1:length(x)
      y(iDim,i) = mean(yrep(whichBin==i));
    end
  end
end
end

function [nData, nSamples] = getNormalizedBinnedData2D(data, x)
for iDim=1:2
  if isfield(data, 'errors')
    nData(iDim,:) = hist(data.errors(iDim,:), x(:,iDim))';
    nData(iDim,:) = nData(iDim,:) ./ sum(nData(iDim,:),'all');
    nSamples = size(data.errors,2);
  else
    for i=1:length(x)
      distM(:,i) = (data.changeSize(iDim,:) - x(i,iDim)).^2;
    end
    [tmp, whichBin] = min(distM,[],2);
    for i=1:length(x)
      nData(iDim,i) = mean(data.afcCorrect(whichBin==i));
    end
    nSamples = numel(data.changeSize(iDim,:));
  end
end
end
