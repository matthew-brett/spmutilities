function [bvals, bvecs] = cbu_dti_params(input_files, is_spm_converted)
% Find bvals and bvecs from DICOM file list or already read headers
% Using Guy William's algorithm
%
% input_files   - string specifying dicom DTI files or
%                     pre-read dicom headers (from spm_dicom_headers)
%                 If not specified, fetches from GUI
% 
% bvals   - N x 1 vector of B values for each of the N dicom DTI files
% bvecs   - N x 3 matrix, with X Y Z (voxel coordinates) normalized
%              diffusion direction for each dicom DTI file
% 
% Examples
% Fetch files from GUI return values
% >> [bvals bvecs] = cbu_dti_params()
% Pass string matrix specifying files
% >> [bvals bvecs] = cbu_dti_params(file_list)
%
% Algorithm GW, mistakes MB

if nargin < 1
  input_files = spm_select(Inf, '^.*\.dcm$', 'Select DICOM files');
end
if nargin < 2
  is_spm_converted = 1;
end
% If already read headers, use these
if iscell(input_files) && isstruct(input_files{1})
    hdrs = input_files;
else % read them
  hdrs = spm_dicom_headers(input_files);
end

% Model SPM's conversion flip in Y
if is_spm_converted
  y_flipper = diag([1 -1 1]);
else
  y_flipper = eye(3);
end

% Get voxel to dicom rotation matrix
orient           = reshape(hdrs{1}.ImageOrientationPatient,[3 2]);
orient(:,3)      = null(orient');
if det(orient)<0, orient(:,3) = -orient(:,3); end;
vox_to_dicom = orient;
vox_to_dicom = vox_to_dicom * y_flipper;
mat = inv(vox_to_dicom);

n_hdrs = length(hdrs);
bvals = zeros(n_hdrs, 1);
bvecs = zeros(n_hdrs, 3);
for h = 1:n_hdrs
  % Read B_matrix info
  bm_info = hdrs{h}.CSAImageHeaderInfo(78);
  % If no B_matrix, this is 0 B value
  if length(bm_info.item) == 0
    continue
  end
  bm = str2num([bm_info.item.val]);
  B_matrix = [bm(1:3); bm(2) bm(4) bm(5); bm(3) bm(5) bm(6)];
  % find max eigenvalue, eigenvector from B_matrix
  [vecs vals] = eig(B_matrix);
  vals = max(vals);
  [bvals(h) i] = max(vals);
  dbvec = vecs(:,i);
  % For convenience, turn vectors to point towards positive X
  if dbvec(1) < 0
    dbvec = dbvec * -1;
  end
  bvecs(h,:) = [mat * dbvec]';
end
