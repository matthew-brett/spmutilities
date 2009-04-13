function [F, p, E] = contrast_comparison(SPM, con, modeled_cols, contrast_cols)
% Analyzes SPM file for relationship of contrast and movement parameters
% FORMAT [F, p, E] = contrast_comparison(SPM, modeled_cols, contrast_cols)
%
% The SPM file should be first-level (fMRI), and be already estimated
% This because we need the filtered design matrix
%
% The function requires marsbar to be on the path
%
% Input
% SPM          - SPM.mat file path or structure
% con          - contrast number from SPM.xCon, or vector, to model
% modeled_cols - columns in SPM design to use to model con effect
% contrast_cols - columns in modeled cols to use as contrast
%                 (default is to use all modeled columns)
% Returns
% F            - F statistic for effect in contrast columns
% p            - p value for F statistic

marsbar on
F = [];
p = [];
if ischar(SPM)
  load(SPM)
end
C1 = sf_get_contrast(con, SPM, 'modeled contrast');
if size(C1, 2) > 1
  error('Single row model contrast only, sorry about that')
end
if nargin < 4
  contrast_cols = 1:length(modeled_cols);
end
if isempty(modeled_cols) | isempty(contrast_cols)
  error('Need non-empty modeled and contrast columns');
end

% Get the filtered design matrix
fX = SPM.xX.xKXs.X;

% Multiply by the modeled contrast to get modeled contrast time course
Y = fX * C1;

% Get the model from modeled columns
X = fX(:, modeled_cols);
names = SPM.xX.name(modeled_cols);
% Add (unfiltered) block means
D = mardo(SPM);
br = block_rows(D);
n_blocks = length(br);
Xm = zeros(size(X,1), n_blocks);
for bno = 1:n_blocks
  Xm(br{bno},bno) = 1;
end
X = [X Xm];
names = [names repmat({'block'},1,n_blocks)];

% Make the new F contrast from new model columns
n_cols = length(contrast_cols);
con_con = zeros(size(X, 2), n_cols);
for ccol = 1:n_cols
  con_con(contrast_cols(ccol), ccol) = 1;
end
% Run the stats by chucking into marsbar
mSPM.xX.X = X;
mSPM.xX.name = names;
D = mardo_2(mSPM);
E = estimate(D, Y);
[E Ic] = add_contrasts(E, {'contrast contrast'}, {'F'}, {con_con});
S = compute_contrasts(E, Ic);
F = S.stat;
p = S.P;
return

function C = sf_get_contrast(con_no, SPM, con_type)
% Process contrast inputs
if isempty(con_no)
  [con_no, xCon] = spm_conman(SPM, 'T|F', 1, ['Select ' con_type]);
end
if numel(con_no) == 1 % Contrast number (probably)
  C = xCon(con_no).c;
  return
end
% More than one element
C = con_no;
if size(C, 1) == 1
  C = C';
end
design_colno =  size(SPM.xX.X, 2);
if size(C, 1) ~= design_colno
  msg = fprintf(...
      'Contrast should have same number of rows as design (=%d)', ...
      design_colno)
  error(msg)
end
return


