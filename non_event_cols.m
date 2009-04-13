function [non_ev_cols, eno_cols] = non_event_cols(SPM)
% Return column numbers for non-event columns
% FORMAT [non_ev_cols, eno_cols] = non_event_cols(SPM)
% 
% Input
% SPM       - SPM design structure
% 
% Output
% non_ev_cols    - column numbers for design effects that are not 
%                  event effects, and are covariates of interest
% eno_cols     - columns that model effect of no interest
%
% Requires marsbar on the path

try
  marsbar on
catch
  error('You need marsbar on the path')
end

con = [];
if nargin < 1
  SPM = spm_select([0 1], '^SPM.*\.mat$', 'Select SPM.mat file to analyze');
end
if isempty(SPM)
  return
end

D = mardo(SPM);
if ~is_fmri(D)
  error('Marsbar does not recognize this as an FMRI design');
end

bcs = block_cols(D);
bcs = cat(2, bcs{:});
X = design_matrix(D);
XM = size(X, 2);
ev_col_list = sf_session_event_cols(D);
other_cols = 1:XM;
other_cols(ismember(other_cols, [ev_col_list{:}])) = [];
cols_in_blocks = ismember(other_cols, bcs);
non_ev_cols = other_cols(cols_in_blocks);
eno_cols = other_cols(~cols_in_blocks);
return

function col_list = sf_session_event_cols(D)
% Return the event modeling columns for each session
all_es = event_specs(D);
session_numbers = all_es(1,:);
if any(diff([0 session_numbers])>1)
  error('Out of order sessions');
end
for ss = unique(session_numbers)
  es = all_es(:,all_es(1,:)==ss);
  cols = [];
  for esno = 1:size(es,2)
    cols = [cols event_cols(D, es(:,esno))];
  end
  if ~all(cols == sort(unique(cols)))
    error('Non-unique or out of order columns');
  end
  col_list{ss} = cols;
end
return


