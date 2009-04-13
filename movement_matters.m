function [F, p] = movement_matters(SPM, conno)
% Test whether movement explains contrast of interest
% FORMAT [F, p] = movement_matters(SPM, conno)
%
% Inputs
% SPM            - SPM design structure or file name
% conno          - contrast that might be explained by movement
%
% Outputs
% F              - F value relating contrast to movement
% p              - p value for F value
  
F = [];
p = [];
if nargin < 1
  SPM = spm_select([0 1], '^SPM.*\.mat$', 'Select SPM.mat file to analyze');
end
if isempty(SPM)
  return
end
if nargin < 2
  conno = [];
end

% This tries to guess which of your design columns are 
% movement columns.  In this case, it finds covariates
% of interest (and excludes events)
[movcols blkcols] = non_event_cols(SPM);
if isempty(movcols)
  error('Design does not seem to have movement parameters')
end
con_cols = [movcols blkcols];
con_con = 1:length(movcols);
[F, p] = contrast_comparison(SPM, conno, con_cols, con_con);

