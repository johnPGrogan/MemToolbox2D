% MEMMODELS2D
%
% One component models:
%   AllGuessingModel2D        - only guessing
%   NoGuessingModel2D         - just precision, no guessing
%
% Mixture models for a single set size:
%   StandardMixtureModel2D    - guess rate and precision
%   SwapModel2D               - guess rate, precision and swaps to other items.
%   SwapModelCovar2D          - guess rate, precision (variance + covariance) and swaps to other items.
%   SwapModelSepMisbind2D     - guess rate, precision and swaps to two sets of other items
%   SwapModelSepPrec2D        - guess rate, 2 precisions (e.g. non-probed targets & distractors) and swaps to other items.
%   SwapModelXYVar2D          - guess rate, 2 precisions (X and Y dimensions) and swaps to other items.
%   SwapModelXYCovar2D        - guess rate, 2 precisions (X and Y dimensions) + covariance of precision, and swaps to other items.
%   VariablePrecisionModel2D  - guess rate and variable precision
%   EnsembleIntegrationModel2D - integration with distractors shifts reports
%
% Models parameterized based on set size:
%   SlotModel2D               - capacity and precision (no benefit when cap.>setsize)
%   SlotsPlusResourcesModel2D - capacity and precision (more juice when cap.>setsize)
%   SlotsPlusAveragingModel2D - capacity and precision (more slots/item when cap.>setsize)
%   ContinuousResourceModel2D - capacity juice split among all items equally
%
% Models that depend on delay duration:
%   ExponentialDecayModel2D   - capacity K and sd, plus objects drop out over time
%
% Model wrappers:
%   WithBias2D               - adds two bias terms (mu for X and Y coords) to any model
%   WithRadialBias2D         - adds a bias terms (mu) to any model - radial
%                               bias towards centre/point
%   TwoAFC2D                  - converts a model so that can be fit to 2afc data
%   WithLapses2D              - adds inattentional parameter to any model
%   WithResponseSampling2D    - sample guesses from a response density
%                               function built over all responses from a
%                               person