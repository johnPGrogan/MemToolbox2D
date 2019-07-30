% Install the MemToolbox2D in the current location. Note that this does not
% copy any files -- you must first put MemToolbox where you want it to
% be located, then run this file to add the files to the path.
% 
% This 2D version is adapted from code in the MemToolbox package(Suchow, J. W.,
% Brady, T. F., Fougnie, D., & Alvarez, G. A. (2013). Modeling visual working 
% memory with the MemToolbox. Journal of Vision, 13(10):9, 1–8. 
% doi:10.1167/13.10.9. MemToolbox.org). 
% It was adapted by John P. Grogan, 2019.
%

addpath(genpath(pwd));
rmpath(genpath(fullfile(pwd, '.git')));
if savepath() == 0
  fprintf(['\nMemToolbox2D successfully added to PATH!\n' ...
    'You are ready to go!\n\n' ...
    'Try calling MemFit2D() to get started.\n\n']);
else
  fprintf(['\nMemToolbox2D failed to add itself to the path!\n' ...
    'Probably this means you don''t have write permissions\n' ...
    'for the pathdef.m file. Here is MATLAB''s error:\n\n']);
  savepath();
end
