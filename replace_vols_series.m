function vols = replace_vols_series(vols, bad_vol_numbers, orig_sdir)
% Replace bad volumes with nearest volumes 
% FORMAT vols = replace_vols_series(vols, bad_vol_numbers, orig_sdir)
% 
% Inputs
% vols             - vol structs or image filenames
% bad_vol_numbers  - vector of volume indices 1..N to be replaced
%                    (where N is the number of volumes in vols)
% orig_sdir        - directory to copy bad volumes into as backup
%
% Returns
% vols             - new vol struct time series

if ischar(vols)
  vols = spm_vol(vols);
end
if isempty(bad_vol_numbers)
  warning('No spikes passed')
  return
end
n_vol = numel(vols);
if nargin < 4
  orig_sdir = [];
end
if isempty(orig_sdir)
  orig_sdir = 'discarded_volumes';
end
img_path = fileparts(vols(1).fname);
orig_path = fullfile(img_path, orig_sdir);
if ~exist(orig_path, 'dir')
  mkdir(img_path, orig_sdir);
end
clusts = sf_clusters(bad_vol_numbers);
for cno = 1:numel(clusts)
  C = clusts{cno};
  rep_nos = sf_replace_nos(C, n_vol);
  if isempty(rep_nos)
    error('Cannot find replacement numbers for cluster');
  end
  % To allow some caching of read replacement volumes
  last_rep_no = 0;
  for vno_cno = 1:numel(C)
    vno = C(vno_cno);
    % Copy volume to subdirectory
    vol = vols(vno);
    img = spm_read_vols(vol);
    [pth fn ext] = fileparts(vol.fname);
    nvol = vol;
    nvol.fname = fullfile(orig_path, [fn ext]);
    spm_write_vol(nvol, img);
    % Replace with replacement 
    rep_no = rep_nos(vno_cno);
    rep_vol = vols(rep_no);
    if rep_no ~= last_rep_no
      rep_img = spm_read_vols(rep_vol);
    end
    rep_vol.fname = vol.fname;
    vol = spm_write_vol(rep_vol, rep_img);
    last_rep_no = rep_no;
  end
end
return

function clusts = sf_clusters(bad_nos)
% Breaks series of indices into contiguous sets
% Return these sets as cell array
bad_nos = sort(unique(bad_nos));
clusts = {};
last_bad = -1;
for i = 1:numel(bad_nos)
  bad_no = bad_nos(i);
  if last_bad == bad_no-1
    % continue cluster
    clust = [clust bad_no];
  else
    % new cluster
    if i > 1
      clusts{end+1} = clust;
    end
    clust = [bad_no];
  end
  last_bad = bad_no;
end
clusts{end+1} = clust;
return

function rep_nos = sf_replace_nos(clust, n_vol)
% Takes cluster indices and returns new replacement indices
% If possible the lowest N/2 indices replaced by the index
% to the left, and highest N/2 replaced by that to right
% If cluster indices start at 1, all replaced by right, similarly
% for cluster indices ending at n_vol (number in series)
% Empty matrix returned if cluster is same as 1:n_vol
rep_left = clust(1)-1;
rep_right = clust(end)+1;
if rep_left < 1 
  if rep_right > n_vol
    rep_nos = [];
    return
  end
  rep_left = rep_right;
end
if rep_right > n_vol
  rep_right = rep_left;
end
n_reps = numel(clust);
rep_nos = zeros(1,n_reps) + rep_right;
L = ceil(n_reps / 2);
rep_nos(1:L) = rep_left;
return
