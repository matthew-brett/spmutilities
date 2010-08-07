function varargout = undistort_reslice(phasemap, imgs, options)
% undistorts images based on phasemap
% FORMAT varargout = undistort_reslice(phasemap, imgs, options)
%
% phasemap     - image with unwrapped phase values
% imgs         - images to undistort (in same [mm/matfile] space as phasemap)
% options      - structure with fields:
%                distinfo - specifies the units of the phasemap
%                   'mm' - map is already adjusted to mm distortion for
%                          this EPI sequence
%                   or: 
%                   structure with fields:
%                     evoltime - evolution time in seconds
%                     epi_bandwidth - bandwidth of _this_ EPI
%                interp - interpolation method - see spm_sample_vol
%                outimgs - one more more of:
%                      'r'  resampled
%                      's'  resampled+smoothed
%                      'm'  mean image
%                      (e.g. 'rsm' for all output images)
%                FWHM - fwhm for smoothing if required
%                wrap - 3 vector specifying wrapping for reslice
%                prefix - prefix for output images [u_r] or [wu_r] (norm)
%                normmat - filename or struct with normalization params
%                        (does realignment only if empty, empty is
%                        default)
%                --- following options are specific for normalization ---
%                bb - 2 x 3 matrix specifying required bounding box in mm
%                     (if any not finite, uses template dims) 
%                vox - 1 x 3 vector of voxel sizes in mm
%                     (if any not finite, uses template voxel sizes)
%                preserve - preserve concentration of tissue 
%
% or (barely hidden feature)
% FORMAT opts = undistort_reslice2('def_opts')
% FORMAT opts = undistort_reslice2('def_opts', 'normalize')
% returns the default options structures for you to play with
  
if nargin < 1
  phasemap = '';
end
if strcmp(phasemap, 'def_opts')
  % you asked for default options
  normf = 0;
  if nargin > 1
    if strcmp(imgs, 'normalize')
      normf = 1;
    end
  end
  varargout{1} = get_def_opts(normf);
  return
end
if nargin < 2
  imgs = '';
end
if nargin < 3
  options = [];
end

% Normalizing or not?
normf = isfield(options, 'normmat') & ~isempty(options.normmat);

% Options
options = fillmerge(options, get_def_opts(normf));

% Check input files
if isempty(phasemap) & ~options.no_phasemap
  phasemap = spm_get(1, 'img', 'Select phase map in EPI space');
end
if isempty(imgs)
  imgs = spm_get(Inf, 'img', 'Select images to undistort');
end

% to vols as necessary
if ischar(phasemap) & ~options.no_phasemap
  phasemap = spm_vol(phasemap);
end
if ischar(imgs)
  imgs = spm_vol(imgs);
end

% output filename and descrips
prefix = options.prefix;
descrip = 'undist reslice';

% normalize or no
if normf
  prm = get_prm(options.normmat);
  [def_options.bb, def_options.vox] = bbvox_from_V(...
      struct('mat', prm.MG, 'dim', prm.Dims(1,:)));
  if ~all(finite(options.vox(:))), options.vox = def_options.vox; end;
  if ~all(finite(options.bb(:))),  options.bb  = def_options.bb;  end;
end

% masking question if required
if options.mask == -1 
  options.mask = (spm_input('Mask the resliced images?','+1','y/n') == 'y');
end;

% resolve output flags
smoothf =  any(options.outimgs == 's');
resampf =  any(options.outimgs == 'r');
meanf = any(options.outimgs == 'm');

nimgs = prod(size(imgs));

% output image stuff
tempv = imgs(1);
if normf % normalizing on the way
  % template image from bounding box etc
  tvox = options.vox;
  bb = options.bb;
  tdim = floor(diff(bb)./tvox)+1;
  origin         = round(-bb(1,:)./tvox + 1);
  off            = -tvox.*origin;
  tmat      = [tvox(1) 0 0 off(1); ...
	       0 tvox(2) 0 off(2); ...
	       0 0 tvox(3) off(3); ...
	       0 0 0 1];
  tempv.dim(1:3) = tdim;
  tempv.mat = tmat;
  tempv.descrip = sprintf('norm %s)',descrip);  
else
  % first image as template
  tempv.descrip = descrip;
  tdim = tempv.dim(1:3);
  tmat = tempv.mat;
  tvox = sqrt(sum(tmat(1:3,1:3).^2));
end

% check dim for first image
dim1 = imgs(1).dim(1:3);
mat1 = imgs(1).mat;
vox1 = sqrt(sum(mat1(1:3,1:3).^2));

linfun = inline('fprintf(''%-60s%s'', x,sprintf(''\b'')*ones(1,60))');

nvox = prod(tdim);

% if normalizing, to mm space of first image
if ~normf % not normalizing
  [X Y Z] = ndgrid(1:tdim(1), 1:tdim(2), 1:tdim(3));
  XYZ = tmat * [X(:) Y(:) Z(:) ones(nvox, 1)]';
else % sadly, we _are_ normalizing, and now I'm confused

  if (prod(size(prm.Tr)) == 0) | options.linear_only
      [X Y Z] = ndgrid(1:tdim(1), 1:tdim(2), 1:tdim(3));
      XYZ = [X(:) Y(:) Z(:) ones(nvox, 1)]';
  else % headache - nonlinear    
    linfun('Getting deformation field');

    Vox = tvox;
    Dims = prm.Dims;
    
    % voxels in normalized image, relative to origin
    x = (bb(1,1):Vox(1):bb(2,1))/Dims(3,1) + Dims(4,1);
    y = (bb(1,2):Vox(2):bb(2,2))/Dims(3,2) + Dims(4,2);
    z = (bb(1,3):Vox(3):bb(2,3))/Dims(3,3) + Dims(4,3);
    Dim = [length(x) length(y) length(z)];
    
    % images to store indices
    [X1 Y1 Z1] = deal(zeros(Dim));

    Tr = prm.Tr;
    
    X = x'*ones(1,Dim(2));
    Y = ones(Dim(1),1)*y;
    basX = spm_dctmtx(Dims(1,1),size(Tr,1),x-1);
    basY = spm_dctmtx(Dims(1,2),size(Tr,2),y-1);
    basZ = spm_dctmtx(Dims(1,3),size(Tr,3),z-1);
    for j=1:length(z),   % Cycle over planes
      % --------------
      tx = reshape(reshape(Tr(:,:,:,1),size(Tr,1)*size(Tr,2),size(Tr,3)) *basZ(j,:)', size(Tr,1), size(Tr,2) );
      ty = reshape(reshape(Tr(:,:,:,2),size(Tr,1)*size(Tr,2),size(Tr,3)) *basZ(j,:)', size(Tr,1), size(Tr,2) );
      tz = reshape(reshape(Tr(:,:,:,3),size(Tr,1)*size(Tr,2),size(Tr,3)) *basZ(j,:)', size(Tr,1), size(Tr,2) );
      X1(:,:,j) = X    + basX*tx*basY';
      Y1(:,:,j) = Y    + basX*ty*basY';
      Z1(:,:,j) = z(j) + basX*tz*basY';
    end
    % reassemble to XYZ matrix
    XYZ = [X1(:) Y1(:) Z1(:) ones(nvox,1)]';
    clear X1 Y1 Z1 x y z basX basY basZ tx ty tz
  end
  
  % apply linear transform to mm space of image that was normalized
  XYZ = prm.MF * prm.Affine * XYZ;
end
clear X Y Z

if ~options.no_phasemap
  % resample from phasemap
  iXYZ = phasemap.mat \ XYZ;
  pvals = spm_sample_vol(phasemap, iXYZ(1,:), iXYZ(2,:), iXYZ(3,:), 1);
  
  % phase encode direction
  phed = options.phase_encode_dim;
  
  % adjust according to distortion info
  if ischar(options.distinfo)
    if strcmp(options.distinfo, 'mm')
      % already in mm distortions -> voxels in first image
      pvals = pvals / vox1(phed);
    else
      error(['Unexpected string ' options.distinfo ' for distinfo'])
    end
  elseif isstruct(options.distinfo)
    % unhandled; something to do with evolution time and bandwidth
    % Go to mm from here
  else
    error('Neither string nor struct for distinfo');
  end

  % undistort before realignment
  % to space of first image, apply undistortion, back to mm
  XYZ = mat1 \ XYZ;
  XYZ(phed, :) = XYZ(phed, :) + pvals;
  XYZ = mat1 * XYZ;
  clear pvals
end % if using phasemap

X = XYZ(1,:);
Y = XYZ(2,:);
Z = XYZ(3,:);
clear XYZ

% masking
if options.mask | meanf,
  linfun('Computing mask..');
  spm_progress_bar('Init',nimgs,'Computing available voxels',...
		   'images completed');
  Count    = zeros(tdim);
  if meanf
    meanimg  = zeros(tdim);
  end
  tiny = 5e-2; % From spm_vol_utils.c
  lolim = ones(3,1) - tiny;
  hilim = dim1' + tiny;
  tmp = zeros(1,nvox);
  for i = 1:nimgs
    M = inv(imgs(i).mat);
    X1=M(1,1)*X+M(1,2)*Y+(M(1,3)*Z+M(1,4));
    Y1=M(2,1)*X+M(2,2)*Y+(M(2,3)*Z+M(2,4));
    Z1=M(3,1)*X+M(3,2)*Y+(M(3,3)*Z+M(3,4));

    tmp = tmp + (...
	  X1 > lolim(1) & X1 < hilim(1) & ...
	  Y1 > lolim(2) & Y1 < hilim(2) & ...
	  Z1 > lolim(3) & Z1 < hilim(3)...
	);
    spm_progress_bar('Set',i);
  end;
  Count = reshape(tmp, tdim);
  if options.mask
    msk = find(Count<nimgs);
  end
  clear tmp
end;

% parameters for smoothing
if smoothf
  [sx, sy, sz, ijk] = convparams(options.FWHM, tvox);
  simg = zeros(tdim);
end

spm_progress_bar('Init',length(imgs),'Undistorting Images');
% resample vols
for i = 1:nimgs
  img = imgs(i);
  outv = tempv; 
  if ~all(dim1 == img.dim(1:3))
    error('Images must all have same voxel dimensions');
  end
  M = inv(img.mat);
  X1=M(1,1)*X+M(1,2)*Y+(M(1,3)*Z+M(1,4));
  Y1=M(2,1)*X+M(2,2)*Y+(M(2,3)*Z+M(2,4));
  Z1=M(3,1)*X+M(3,2)*Y+(M(3,3)*Z+M(3,4));
  ivals = spm_sample_vol(img,X1,Y1,Z1,options.interp);
  ivals = reshape(ivals, tdim);
  if meanf
    meanimg = meanimg + ivals;
  end

  [p f e] = fileparts(img.fname);

  if resampf
    if options.mask
      ivals(msk) = NaN;
    end
    outv.fname = fullfile(p, [prefix f e]);
    spm_write_vol(outv, ivals);
  end
  
  if smoothf
    if options.mask
      ivals(msk) = 0;
    end
    spm_conv_vol(ivals, simg, sx, sy, sz, ijk); 
    outv.fname = fullfile(p, ['s' prefix f e]);
    outv.descrip = sprintf('%s - conv (%g,%g,%g)',outv.descrip, options.FWHM);
    spm_write_vol(outv, simg);
  end
  
  spm_progress_bar('Set',i);
end
clear ivals msk simg iXYZ XYZ

if meanf
  meanv = tempv;
  [p f e] = fileparts(meanv.fname);
  meanv.fname = fullfile(p, [prefix 'mean' f e]);
  meanv.descrip = 'spm  - umean image';
  msk = Count > 0;
  meanimg(msk) = meanimg(msk) ./ Count(msk);
  meanv = spm_write_vol(meanv, meanimg);
end

spm_progress_bar('Clear');
return

% Subfunctions 

%_______________________________________________________________________
function [Dims,Affine,MF,MG,Tr] = load_params(matname)
load(deblank(matname))
if (exist('mgc') ~= 1)
	error(['Matrix file ' matname ' is the wrong type.']);
end
if (mgc ~= 960209)
	error(['Matrix file ' matname ' is the wrong type.']);
end

Tr = reshape(Transform,[Dims(2,:) 3]);
return;
%_______________________________________________________________________


function s1 = fillmerge(s1, s2)
% fills missing or empty fields from s1 with nonmissing fields in s2
fns = fieldnames(s2);
for f = 1:length(fns)
  fn = fns{f};
  if ~isfield(s1, fn) | isempty(getfield(s1, fn))
    s1 = setfield(s1, fn, getfield(s2, fn));
  end
end
return

function [x, y, z, ijk] = convparams(FWHM, VOX)
% compute parameters for spm_conv_vol
%-----------------------------------------------------------------------
s  = FWHM./VOX;					% voxel anisotropy
s  = max(s,ones(size(s)));			% lower bound on FWHM
s  = s/sqrt(8*log(2));				% FWHM -> Gaussian parameter

x  = round(6*s(1)); x = [-x:x];
y  = round(6*s(2)); y = [-y:y];
z  = round(6*s(3)); z = [-z:z];
x  = exp(-(x).^2/(2*(s(1)).^2)); 
y  = exp(-(y).^2/(2*(s(2)).^2)); 
z  = exp(-(z).^2/(2*(s(3)).^2));
x  = x/sum(x);
y  = y/sum(y);
z  = z/sum(z);

ijk  = -([length(x) length(y) length(z)]-1) /2  ;
return

%_______________________________________________________________________
function [bb,vx] = bbvox_from_V(V)
vx = sqrt(sum(V.mat(1:3,1:3).^2));
if det(V.mat(1:3,1:3))<0, vx(1) = -vx(1); end;

o  = V.mat\[0 0 0 1]';
o  = o(1:3)';
bb = [-vx.*(o-1) ; vx.*(V.dim(1:3)-o)]; 
return;
%_______________________________________________________________________

function opts = get_def_opts(normf)
% get default options

opts = struct(...
    'distinfo', 'mm',... % phasemap in terms of mm distortion
    'interp', 1, ...     % - resampling value (see spm_sample_vol.m) [1]
    'phase_encode_dim', 2,... % phase encode dimension [2] (= Y dimension)
    'mask', 1,...        % - 0 or 1, mask images a la spm_reslice [1]
    'outimgs', 'r',...   % output images
    'FWHM', 8,...        % smoothing FWHM for smoothed images
    'prefix', 'u_r',...  % prefix for output images
    'normmat', '',...    % normalization structure or matrix
    'linear_only', 0,... % flag to use linear norm component only 
    'no_phasemap', 0 ... % don't use phasemap (for testing)
    );

% get SPM99 defaults
global sptl_MskOptn sptl_BB sptl_Vx
all_defs = struct(...
    'mask', sptl_MskOptn,...
    'bb',   sptl_BB,...
    'vox',  sptl_Vx);
if normf
  opts.prefix = ['n' opts.prefix];
end
opts = fillmerge(all_defs, opts);
return

function prm = get_prm(normmat)
if ischar(normmat), prm = load(normmat); end
% Try to import SPM2 normalization
if ~isfield(prm, 'mgc')
  tmp = prm;
  prm.MG = tmp.VG.mat;
  prm.MF = tmp.VF.mat;
  prm.Affine = tmp.Affine;
  prm.Tr = tmp.Tr;
  sz = [size(prm.Tr) 0 0 0];
  vox  = sqrt(sum(prm.MG(1:3,1:3).^2));
  origin = prm.MG\[0 0 0 1]';
  origin = round(origin(1:3)');
  prm.Dims = [tmp.VG.dim(1:3);...
	      max([0 0 0; sz(1:3)]);...
	      vox;...
	      origin;...
	      tmp.VF.dim(1:3);...
	      sqrt(sum(prm.MF(1:3,1:3).^2))];
else
  prm.Tr = reshape(prm.Transform,[prm.Dims(2,:) 3]);
end
return
