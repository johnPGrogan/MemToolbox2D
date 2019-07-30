function dimensions = GetModelDims(data)
% dimensions = getModelDims(data)
% return data.dimensions if present, else use default of 1366*768


if isfield(data,'dimensions') % use dimensions if supplied, else assume screen is [1366 * 768]
  dimensions = [data.dimensions(1); data.dimensions(2)];
else
  dimensions = [1366; 768];
end
  
end