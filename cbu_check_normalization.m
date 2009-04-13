function cbu_check_normalization(P, temp)
% FORMAT checknorm(P,temp)
% SPM99 etc function to check normalisation
% 
% To use, start spm 99, then type checknorm in the matlab window
% Select the normalised images, then select the template to 
% check against. Checknorm displays the files one by one
% alongside the template, using the check reg facility from
% spm99.  To see the next image, click on Yes for continue
% in the SPM input window.  The current image name appears
% at the top of the same window.
%
% Matthew Brett 6/99
  
if nargin < 1
  P = spm_get(Inf, 'nmean*.img', 'Files to check against template');
end
if nargin < 2
  temp = spm_get(1, 'img', 'Template image to check against', ...
		 [spm('Dir') filesep 'templates']);
end

fnames = spm_str_manip(P, 'c');
for i=1:size(P,1)
  filen = deblank(fnames(i,:));
  spm_input(filen,1,'d', 'File;');
  spm_check_registration(str2mat(P(i, :), temp));
  if (i<size(P,1)) 
    if ~spm_input('Continue?',2, 'y/n',[1 0],1)
      break
    end
  end
end

