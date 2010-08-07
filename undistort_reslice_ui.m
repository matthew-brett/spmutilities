function undistort_reslice_ui(pfile)
% little GUI for undistort_realign routine

if nargin > 0
  if ischar(pfile)
    pfile = load(pfile);
  end
  [errf phasemap imgs matname] = fieldwarn(pfile, ...
					   {'phasemap','imgs','matname'});
  if errf
    error(phasemap);
  end
  n = length(imgs);
else
  n = spm_input('number of subjects', 1, 'e', 1);
  if (n < 1)
    spm_figure('Clear','Interactive');
    return;
  end
  filt = '*_sn.mat'; end
  for s = 1:n
    matname{s} = spm_get(1, filt, 'Select normalization params');
    phasemap{s} = spm_get(1, ...
			  '*rpm_*_rot_dis.img', ...
			  'Select distortion phasemap');
    
    imgs{s} = spm_get(Inf, 'img', 'Select images to normalize');
  end
end
nopts = undistort_reslice('def_opts', 'normalize');
outimgs = spm_input('Images to output', '+1', 'm', ...
		    'normalized|smoothed normalized|both', [1 2 3], 3);
if outimgs ~= 2 % normalized
  nopts.outimgs = 'r';
else
  nopts.outimgs = '';
end
if outimgs > 1 % smoothed
  nopts.FWHM = spm_input('FWHM for smoothing', '+1');
  nopts.outimgs = [nopts.outimgs 's'];
end
for s = 1:n
  nopts.normmat = matname{s};
  undistort_reslice(phasemap{s}, imgs{s}, nopts);
end
save undistort_reslice.mat phasemap imgs matname 


function [errf,varargout] = fieldwarn(st, fieldn)
errf = 0;
msg = '';
fno = length(fieldn);
varargout = cell(1, fno);
for f = 1:fno
  fn = fieldn{f};
  if ~isfield(st,fn)
    errf = 1;
    msg = strvcat(msg, ['Expecting field: ' fn]);
  else
    varargout{f} = getfield(st, fn);
  end
end
if errf, varargout{1} = msg; end
return
