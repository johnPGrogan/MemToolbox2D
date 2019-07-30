% DOESMODELREQUIREEXTRAINFO2D checks if a 2D model pdf requires more than 
% data.errors
%
% r = DoesModelRequireExtraInfo2D(model)
%
% This function is a helper function used by functions that attempt to
% evaluate a model's pdf at specific values independent of the specified
% data. It check if a model.pdf function requires more than just specifying
% the data.errors to evaluate at (for example, the SwapModel2D also requires
% you to include data.distractors, and thus will return true).
%
% This 2D version is adapted from code in the MemToolbox package(Suchow, J. W.,
% Brady, T. F., Fougnie, D., & Alvarez, G. A. (2013). Modeling visual working 
% memory with the MemToolbox. Journal of Vision, 13(10):9, 1–8. 
% doi:10.1167/13.10.9. MemToolbox.org). 
% It was adapted by John P. Grogan, 2019.
%

function r = DoesModelRequireExtraInfo2D(model)
  r = false;
  try
    data.errors = [-1 0 1;-1 0 1];
    params = num2cell(model.start(1,:));
    model.pdf(data, params{:});
  catch
    r = true;
  end
end
