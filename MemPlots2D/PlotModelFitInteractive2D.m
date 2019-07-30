% PlotModelFitInteractive2D allow manipulation of the fit with sliders
% This functions plots the model's probability density function
% overlaid on a histogram of the data. The plot is interactive, with
% a slider that allows you to adjust each of the parameters of the model
% and see the impact on the pdf.
%
%    figHand = PlotModelFitInteractive2D(model, params, data, varargin)
%
% 'params' can be either a maxPosterior, a fullPosterior or a
% posteriorSamples.
%
% Optional parameters:
%  'PdfColor' - the color to plot the model fit with. Note that this color
%  is automatically 'faded' if the model fit being shown is not the max
%  posterior fit.
%
%  'MarginalPlots' - If this is set to true, then under each slider bar is
%  plotted the current likelihood function for each parameter, conditioned
%  on the values of the other parameters.
%
%  'NewFigure' - whether to make a new figure or plot into the currently
%  active subplot. Default is false (e.g., plot into current plot).
%
% 
% This 2D version is adapted from code in the MemToolbox package(Suchow, J. W.,
% Brady, T. F., Fougnie, D., & Alvarez, G. A. (2013). Modeling visual working 
% memory with the MemToolbox. Journal of Vision, 13(10):9, 1–8. 
% doi:10.1167/13.10.9. MemToolbox.org). 
% It was adapted by John P. Grogan, 2019.
%

function figHand = PlotModelFitInteractive2D(model, params, data, varargin)
  % Extra parameters
  args = struct('MarginalPlots', false, 'NewFigure', true, ...
    'PdfColor', [0.54, 0.61, 0.06]);
  args = parseargs(varargin, args);
  if args.NewFigure, figHand = figure(); else figHand = []; end

  % If you pass a 'posteriorSamples' struct instead of params
  if isstruct(params) && isfield(params, 'vals')
    params = MCMCSummarize(params, 'maxPosterior');
  end

  % If params has a .logLikeMatrix, assume they passed a fullPosterior from
  % GridSearch
  if isstruct(params) && isfield(params, 'logLikeMatrix')
    posteriorSamples = SampleFromPosterior(params, 500);
    params = MCMCSummarize(params, 'maxPosterior');
  end

  % Ensure there is a model.prior, model.logpdf and model.pdf
  model = EnsureAllModelMethods(model);

  % Initial plot
  paramsCur = params;
  paramsCell = num2cell(params);
  mapLikeVal = model.logpdf(data, paramsCell{:});

  PlotModelFit2D(model, paramsCur, data, 'ShowNumbers', false, ...
        'PdfColor', args.PdfColor,'SubplotDims',[2 2]);
%   pos = get(gca, 'Position');
% 
%   % Decide on spacing of plot, based on how many parameters (and thus
%   % sliders) we need
  Nparams = length(model.paramNames);
  vertSpacing = (.35/Nparams);
  if vertSpacing>0.10
    vertSpacing = 0.10;
  end
  maxPos = min([.35 vertSpacing*Nparams]);
%   pos(2) = maxPos + 0.12;
%   pos(4) = 1 - maxPos - 0.15;
  poss = [.13 .11 .3347 .3412; .5703 .11 .3347 .3412];
  
for iDim = 1:2
    pos = poss(iDim,:);
%   set(gca, 'Position', pos);
    subplot(2,2,iDim+2)
    axis('off')
  height = vertSpacing-0.03;
  mainAxis = gca;

  % If we're also going to plot marginals, make scroll bars shorter
  if args.MarginalPlots
    height=height-0.02;
  end

  % For each parameter, decide how its slider should map to its parameter
  % range, and then make the slider...
  for i=1:Nparams
    if ~isinf(model.upperbound(i))
      MappingFunction{i,iDim} = @(percent) (model.lowerbound(i) ...
        + (model.upperbound(i)-model.lowerbound(i)).*percent);
      InverseMappingFunction{i,iDim} = @(val) ((val-model.lowerbound(i)) ...
        / (model.upperbound(i)-model.lowerbound(i)));
    else
      if isinf(model.lowerbound(i))
        % Should probably use a logistic with variable mean. For now just
        % error!
        error('Can''t have lower and upperbound of a parameter be Inf');
      else
        MappingFunction{i,iDim} = @(percent) (-log(1-percent).*params(i)*2);
        InverseMappingFunction{i,iDim} = @(val) (1-exp(-val/(params(i)*2)));
      end
    end

    invertedMapping = InverseMappingFunction{i,iDim}(params(i));
    slider(i,iDim) = uicontrol(...
      'Parent',gcf,...
      'Units','normalized',...
      'Callback', @(hObject,eventdata) slider_Callback(hObject, i),...
      'Position',[pos(1) maxPos-(i-1)*vertSpacing-height .3 height],...
      'Style','slider', ...
      'UserData', MappingFunction{i}, ...
      'Value', invertedMapping);

    % Create plots for marginals
    if args.MarginalPlots
      myHeight = height+0.03;
      marginalPlot(i,iDim) = axes('Position', ...
        [0.10 maxPos-(i-1)*vertSpacing-myHeight 0.75 myHeight]);
      colormap(palettablecolormap());
      valuesForMarginalPlot{i,iDim} =  MappingFunction{i}(0:0.01:1);
      set(marginalPlot(i,iDim), 'Units', 'pixels');
      newPos = get(marginalPlot(i,iDim), 'Position');
      set(marginalPlot(i,iDim), 'Position', [newPos(1)+15, newPos(2), newPos(3)-30, newPos(4)]);
    end

    uicontrol(...
      'Parent',gcf,...
      'Units','normalized',...
      'BackgroundColor',[1 1 1],...
      'Style','text',...
      'FontSize', 12, ...
      'FontWeight', 'bold', ...
      'Position',[0.02 maxPos-(i-1)*vertSpacing-height 0.13 height],...
      'String', model.paramNames{i});

    curVals(i,iDim) = uicontrol(...
      'Parent',gcf,...
      'Units','normalized',...
      'BackgroundColor',[1 1 1],...
      'Style','edit',...
      'FontSize', 12, ...
      'FontWeight', 'bold', ...
      'Position',[pos(1)+.3 maxPos-(i-1)*vertSpacing-height 0.13 height],...
      'String', sprintf('%0.2f', paramsCur(i)), ...
      'Callback', @(hObject,eventdata) edit_Callback(hObject, i), ...
      'UserData', InverseMappingFunction{i});
  end
end
  % Allow the user to limit this figure to any subset of the data
  if ~isempty(figHand)
    CreateMenus(data, @redrawFig);
  end
  function redrawFig(whichField, whichValue)
    if strcmp(whichField, 'all')
      subplot(1,1,1);
      PlotModelFitInteractive2D(model, params, data, ...
        'MarginalPlots', args.MarginalPlots, 'NewFigure', false, ...
        'PdfColor', args.PdfColor);
    elseif sum(ismember(data.(whichField),whichValue)) > 0
      [datasets,conditionOrder] = SplitDataByField(data, whichField);
      newData = datasets{ismember(conditionOrder,whichValue)};
      subplot(1,1,1);
      PlotModelFitInteractive2D(model, params, newData, ...
        'MarginalPlots', args.MarginalPlots, 'NewFigure', false, ...
        'PdfColor', args.PdfColor);
    end
  end

  PlotMarginals();


  % To draw the conditional distribution of each parameter given the other
  % values of the other parameters (if option is chosen)
  function PlotMarginals()
    if args.MarginalPlots
      for i=1:Nparams
        axes(marginalPlot(i));
        newParamsCell = num2cell(paramsCur);
        for p=1:101 % 0:0.01:1
          newParamsCell{i} = valuesForMarginalPlot{i}(p);
          loglike(p) = model.logpdf(data, newParamsCell{:});
        end
        marginal = exp(loglike-max(loglike));
        imagesc(marginal./nansum(marginal(:)), [0 1]);
        set(gca, 'XTick', [], 'YTick', []);
        set(gca, 'XColor', [.8 .8 .8], 'YColor', [.8 .8 .8]);
      end
    end
  end

  % If they type a new number, set the slider to that value
  function edit_Callback(hObject, which)
    curValue = str2double(get(hObject,'String'));
    if isnan(curValue)
      set(hObject, 'String', sprintf('%0.2f', paramsCur(which)));
      beep; warning('Not a number!');
      return;
    end
    paramsCur(which) = curValue;

    inverseMappingFunc = get(hObject, 'UserData');
    set(slider(which), 'Value', inverseMappingFunc(curValue));
    slider_Callback(slider(which), which);
  end

  % When the slider moves, find out what value for the parameter that
  % should correspond to, and set the edit box and the plot to show that
  function slider_Callback(hObject, which)
      
    % find which subplot clicked on
    axesClicked = get(gcf,'CurrentObject');
    axToUpdate = 2 - (axesClicked.Position(1) < .5);
    
    curValue = get(hObject,'Value');
    mappingFunc = get(hObject, 'UserData');
    paramsCur(which) = mappingFunc(curValue);
    axes(mainAxis); hold off;
    paramsCell = num2cell(paramsCur);
    curLike = model.logpdf(data, paramsCell{:});
    PlotModelFit2D(model, paramsCur, data, 'ShowNumbers', false, 'PdfColor', ...
      fade(args.PdfColor, exp(curLike - mapLikeVal)),'SubplotDims',[2 2],'SubplotInds',axToUpdate);
    PlotMarginals();
    set(curVals(which,axToUpdate), 'String', sprintf('%0.2f', paramsCur(which)));
%     set(gcf,'Cur
  end
end




