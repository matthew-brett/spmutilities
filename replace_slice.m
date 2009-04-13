function V = replace_slice(vol, slice_number, vol_left, vol_right, prefix)
if nargin < 3
  vol_left = [];
end
if nargin < 4
  vol_right = [];
end
if nargin < 5
  prefix = [];
end
if isempty(prefix)
  prefix = 'fixed_slice_';
end
if isempty(vol_left) & isempty(vol_right)
  error('No volumes to replace with');
end
if ischar(vol)
  vol = spm_vol(vol);
end
if ischar(vol_left)
  vol_left = spm_vol(vol_left);
end
if ischar(vol_right)
  vol_right = spm_vol(vol_right);
end
if slice_number < 1 | slice_number > vol.dim(3)
  error('Slice out of range');
end
img = spm_read_vols(vol);
S = zeros(vol.dim(1:2));
n = 0;
if ~isempty(vol_left)
  img2 = spm_read_vols(vol_left);
  S = S + img2(:,:,slice_number);
  n = n+1;
end
if ~isempty(vol_right)
  img2 = spm_read_vols(vol_right);
  S = S + img2(:,:,slice_number);
  n = n+1;
end
S = S / n;
img(:, :, slice_number) = S;
[path fname ext] = fileparts(vol.fname);
V = vol;
V.fname = fullfile(path, [prefix fname ext]);
V = spm_write_vol(V, img);

