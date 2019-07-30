%PLOTDATA2D plots a histogram of the data
% Can plot either continuous report data, or binned bar graph for 2AFC data
%
%  figHand = PlotData2D(data, [optionalParameters])
%
% Optional parameters:
%  'NumberOfBins' - the number of bins to use in display the data. Default
%  40.
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

function figHand = PlotData2D(data, varargin)
  % Extra arguments and parsing
  args = struct('NumberOfBins', 40, 'NewFigure', true);
  args = parseargs(varargin, args);
  if args.NewFigure, figHand = figure(); else figHand = []; end

  % Clean data up if it is just errors
  if(~isfield(data,'errors')) && (~isfield(data,'afcCorrect'))
    data = struct('errors',data);
  end

  dims = GetModelDims(data); % get the model dimensions
  xLabs = {'X', 'Y'}; % dimensions

  if isfield(data, 'errors')
    % Plot data histogram for continuous report data
      for i = 1:2
        x(i,:) = linspace(-dims(i), dims(i), args.NumberOfBins+1)';

        subplot(1,2,i)
        n(i,:) = hist(data.errors(i,:)', x(i,:));
        bar(x(i,:), n(i,:)./sum(n(i,:)), 'EdgeColor', [1 1 1], 'FaceColor', [.8 .8 .8]);

        hold on;
        set(gca, 'box', 'off');
        xlabel([xLabs{i} ' Error'], 'FontSize', 14);
        ylabel('Probability', 'FontSize', 14);
        topOfY = max(n(i,:)./sum(n(i,:)))*1.20;
        ylim([0 topOfY]);
      end
    
  else

    % Plot binned data for 2AFC data
    for i = 1:2
        subplot(1,2,i)
        set(gcf, 'Color', [1 1 1]);
        x = linspace(-dims(i), dims(i), args.NumberOfBins)';
        for j=2:length(x)
            which = all(data.changeSize>=x(j-1) & data.changeSize<x(j));
            mn(j-1) = mean(data.afcCorrect(which));
            se(j-1) = std(data.afcCorrect(which))./sqrt(sum(which));
        end
        binX = (x(1:end-1) + x(2:end))/2;
        bar(binX, mn, 'EdgeColor', [1 1 1], 'FaceColor', [.8 .8 .8]);
        hold on;
        errorbar(binX, mn, se, '.', 'Color', [.5 .5 .5]);
        xlim([-dims(i) dims(i)]);
        set(gca, 'box', 'off');
        xlabel([xLabs{i} ' Distance'], 'FontSize', 14);
        ylabel('Probability Correct', 'FontSize', 14);
        ylim([0 1]);
    end
  end

  % Allow the user to limit this figure to any subset of the data
  if ~isempty(figHand)
    CreateMenus(data, @redrawFig);
  end
  function redrawFig(whichField, whichValue)
    if strcmp(whichField, 'all')
      cla;
      PlotData2D(data, 'NewFigure', false, 'NumberOfBins', args.NumberOfBins);
    elseif sum(ismember(data.(whichField),whichValue)) > 0
      [datasets,conditionOrder] = SplitDataByField(data, whichField);
      newData = datasets{ismember(conditionOrder,whichValue)};
      cla;
      PlotData2D(newData, 'NewFigure', false, 'NumberOfBins', args.NumberOfBins);
    end
  end
end
