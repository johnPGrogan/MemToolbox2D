% NOGUESSINGMODEL2D returns a structure for a single component model.
% This is the same as StandardMixtureModel2D, but without a guess state.
% The probability distribution is a uniform distribution of error.
%
% This 2D version is adapted from code in the MemToolbox package(Suchow, J. W.,
% Brady, T. F., Fougnie, D., & Alvarez, G. A. (2013). Modeling visual working 
% memory with the MemToolbox. Journal of Vision, 13(10):9, 1–8. 
% doi:10.1167/13.10.9. MemToolbox.org). 
% It was adapted by John P. Grogan, 2019.
%

function model = NoGuessingModel2D()
  model.name = 'No guessing model2D';
	model.paramNames = {'sd'};
	model.lowerbound = 0; % Lower bounds for the parameters
	model.upperbound = Inf; % Upper bounds for the parameters
	model.movestd = 0.5;
	model.pdf = @(data, sd) (mvnpdf(data.errors', [0 0],[sd.^2, sd.^2]));
	model.start = [ 3;  % sd
                 15;  % sd
                 75]; % sd

  % To specify a prior probability distribution, change and uncomment
  % the following line, where p is a vector of parameter values, arranged
  % in the same order that they appear in model.paramNames:
  % model.prior = @(p) (1);

end
