function targResps = ResampleOrBound(targResps, targs, varVec, model, dimensions)
% function targResps = ResampleOrBound(targResps, targs, varVec, model, dimensions)
%
% check whether you wish to resample the generative model or simply set
% items occuring outside the screen to the edges
% if the former, it resamples using mvnrnd up to maxIter times, and gives
% and error if it exceeds that
%
% Inputs: 	targResps are the [x;y] coordinates of responses
%			targs are the [x;y] coordinates of target locations
% 			varVec is a vector of variances (assuming zero covariance) for the mvnrnd()
% 			model is the model used, which can contain 'boundary' and 'maxIter' fields
% 				boundary=1 means items outside screen are moved to closest point on boundary
%					else they are resampled using mvnrnd
% 				maxIter is the max number of iterations for resampling (default=100)
% 			dimensions is the [x;y] dimensions of the screen
%
% Output:   targResps are the new [x;y] coordiantes of sampled responses
% 
% Written by John P Grogan, 2019.


  n = size(targs, 2); % number of trials
  
  toReplace = targResps < 0 | targResps > dimensions'; % find responses outside of screen
  
  if isfield(model, 'boundary') && model.boundary % if edge-constraining chosen
      
      for iDim = 1:2 % for X and Y, move to boundary
        targResps(targResps(:,iDim) > dimensions(iDim),iDim) = dimensions(iDim); 
        targResps(targResps(:,iDim) <0,iDim) = 0;
      end
      
  else % if edge-resampling
      
      % if no maxIterations given, set as 100
      if ~isfield(model,'maxIter'), maxIter=100; else maxIter=model.maxIter; end
      
      nIter=1;
      toReplace = any(toReplace,2); % which trials need to be resampled?
      
      while any(toReplace) % keep going until all are within screen
          targResps1 = mvnrnd(targs', varVec, n); % redraw sample
          targResps(toReplace,:) = targResps1(toReplace,:); % replace
          toReplace = any(targResps < 0 | targResps > dimensions',2); % check which are still outside, repeat if necessary
          nIter=nIter+1; 
          if nIter > maxIter
              error('max iterations for resampling exceeded')
          end
      end
      
  end
  
end