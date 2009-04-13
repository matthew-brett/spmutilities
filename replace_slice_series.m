function vols = replace_slice_series(vols, bad_vol_number, slice_number, prefix)
if ischar(vols)
  vols = spm_vols(vols);
end
if nargin < 4
  prefix = [];
end
vol = vols(bad_vol);
if bad_vol > 1
  left_vol = vols(bad_vol-1);
else
  left_vol = [];
end
if bad_vol < length(vols)
  right_vol = vols(bad_vol+1);
else
  right_vol = [];
end
V2 = replace_slice(vol, slice_number, left_vol, right_vol, prefix);
vols(bad_vol) = V2;
