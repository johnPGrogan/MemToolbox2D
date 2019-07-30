% VARIABLEPRECISIONMODEL2D returns a structure for an infinite scale mixture.
% In such a model, the standard deviations of observers' reports are assumed
% to be themselves drawn from a higher-order variability distribution,
% rather than always fixed.
%
% The default model assumes observers' standard deviations are distributed
% according to a Gaussian distribution.  However, the function takes an
% optional argument 'HigherOrderDist'. VariablePrecisionModel2D('HigherOrderDist',
% 'GammaSD') returns a model where the higher-order distribution of SD is assumed
% to be Gamma, rather than Gaussian. VariablePrecisionModel2D('HigherOrderDist',
% 'GammaPrecision') assumes a distribution over precisions (1/Variance) that is
% Gamma, as in van den Berg et al. (2012).
%
% This 2D version is adapted from code in the MemToolbox package(Suchow, J. W.,
% Brady, T. F., Fougnie, D., & Alvarez, G. A. (2013). Modeling visual working 
% memory with the MemToolbox. Journal of Vision, 13(10):9, 1–8. 
% doi:10.1167/13.10.9. MemToolbox.org). 
% It was adapted by John P. Grogan, 2019.
%

function model = VariablePrecisionModel2D(varargin)
  % Default: Don't include a bias term, and a Gaussian over SDs as
  % higher-order distribution
  args = struct('HigherOrderDist', 'GaussianSD');
  args = parseargs(varargin, args);

  if strcmp(args.HigherOrderDist, 'GaussianSD')
    model = VariablePrecisionModel2D_GaussianSD();
  elseif strcmp(args.HigherOrderDist, 'GammaSD')
    model = VariablePrecisionModel2D_GammaSD();
  elseif strcmp(args.HigherOrderDist, 'GammaPrecision')
    model = VariablePrecisionModel2D_GammaPrecision();
  end
end

