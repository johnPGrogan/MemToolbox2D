function data = ChangeMTB2DUnits(data, denominator, fieldsToConvert)
% change all coord units in a MemToolbox2D data structure
% Can use data.dimensions to normalise to 0-1 (or data.dimensions/100 for
% 0-100), or can use pixelsPerMm to convert into mm (and thereafter into 
% visual degrees.
% 
% Inputs:
%   data = MemToolbox2D data structure
%   denominator = factor to divide by to change units. either [x; y] or
%     scalar. (to make proportion of dimensions use data.dimensions, to
%     make percentage of dimensions, use data.dimesions/100).
%   fieldsToConvert = cell array of fieldnames to convert. default is all
%       but data.n and data.whichIsTestItem
% Outputs:
%   data = data structure in new units
% 
% John Grogan, 2021.
% 

%% check inputs

if ~isstruct(data)
    error('data should be a MemToolbox2D data structure')
end
% get fieldnames
fieldNames = fieldnames(data);

if ~ismatrix(denominator) || numel(denominator) >2
    error('denominator should be scalar or 2x1 column vector')
elseif all(size(denominator) == [1 2])
    denominator = denominator'; % rotate into column
end
if exist('fieldsToConvert','var')
    if ~iscell(fieldsToConvert)
        error('fieldsToConvert must be a cell array of strings')
    end
else
    fieldsToConvert = fieldNames(~ismember( fieldNames,{'whichIsTestItem','n'})); % all but n and whichIsTestItem
end


%% use common names between fields and fieldsToConvert

fieldNames = fieldNames(ismember(fieldNames, fieldsToConvert));

%% convert them
for i = 1:length(fieldNames)
    x = data.(fieldNames{i}); % copy daty
    
    r = size(x,1) ./ size(denominator,1);
    if mod(r, 1) == 0 % if is has even rows
        denom1 = repmat(denominator, r, 1); % repeat rows to match size       
        
        data.(fieldNames{i}) = x ./ denom1;
    
    else
        error('number of rows in field does not allow division by denominator');
    end
    
    
end
