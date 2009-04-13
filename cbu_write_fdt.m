function cbu_write_fdt(bvals, bvecs, outpath, prefix)
% Write FSL FDT bvals and bvecs files from parameters
% FORMAT write_fdt_files(bvals, bvecs, outpath)
% 
% bvals   - Nx1 vector with diffusion gradient strength
%           (one per diffusion acquisition, N=no of acquisitions)
% bvecs   - Nx3 matrix with diffusion gradient directions
%           
% outpath - path to write FDT bvals, bvecs text files
%           (defaults to current working directory)
% prefix  - prefix for bvals, bvecs files in directory
%
% If no input args, fetch bvals, bvecs using cbu_dti_params
% 
% examples
%
% Select files via GUI write bvecs, bvals to current directory
% >> cbu_write_fdt(); 
% Write pre-existing variables to file in given directory
% >> cbu_write_fdt(bvals, bvecs, '/my/dti/analysis/directory');
%
% Matthew Brett 30 April 2007, 17 July 2007
  
if nargin < 2
  bvals = [];
  bvecs = [];
end
if isempty(bvals) | isempty(bvecs)
  [bvals, bvecs] = cbu_dti_params;
end
if nargin < 3
  outpath = pwd;
end
if nargin < 4
  prefix = [];
end
fid = fopen(fullfile(outpath, [prefix 'bvals']), 'wt');
if fid == -1
  error('Could not open bvals output file')
end
fprintf(fid, '   %e', bvals);
fprintf(fid, '\n');
fclose(fid);

fid = fopen(fullfile(outpath, [prefix 'bvecs']), 'wt');
if fid == -1
  error('Could not open bvecs output file')
end
for D = 1:3
  fprintf(fid, '   %e', bvecs(:,D));
  fprintf(fid, '\n');
end
fclose(fid);
return
