% VALIDATEDATA2D checks to make sure that the data is in the expected format
%
% [data, pass] = ValidateData2D(data)
%
% e.g. it checks to see if there are 2 rows of data [x; y] and
% e.g., it checks if it is in the range [-180,180]. if unsalvageable, it
% throws errors. Otherwise, throws warnings and does its best to massage
% data into the range (-180, 180).
%
function [data, pass] = ValidateData2D(data)
    pass = true; % always pass if you make it through without an error()

    % Rename according to MTB standards when appropriate
    if ~isfield(data, 'errors') && isfield(data, 'error')
      data.errors = data.error;
      data = rmfield(data, 'error');
    end

    if(~isDataStruct(data))
      error('Data should be passed in as a struct with a field data.errors or data.afcCorrect');
    end

    if isfield(data, 'errors') && isfield(data, 'afcCorrect')
      error('Your data struct specified both a .errors and a .afcCorrect. Please pass only one kind of data at a time.');
    end
    
    dims = GetModelDims(data);

    % Check that the error values are in the correct range, otherwise massage
    if isfield(data, 'errors')
      if(isempty(data.errors))
        error('The data vector should not be empty.');
      elseif(~isnumeric(data.errors))
        throwRangeError();
      elseif any(data.errors < - dims | data.errors > dims,'all') % errors outside dims
        throwRangeErrors();
      end
    end

end


function throwRangeError()
    error('Data should be in the range (-dims:dims)');
end

% Is the object an MTB data struct? passes iff the object is a struct
% containing a field called 'errors'.
function pass = isDataStruct(object)
  pass = (isstruct(object) && (isfield(object,'errors') || isfield(object, 'afcCorrect')));
end

