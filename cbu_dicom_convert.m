function [filenames, hdrs] = cbu_dicom_convert(dicom_directory, output_directory, filt, output_format)
% Reads DICOM directory, write converted images, return vols, hdrs
% FORMAT [vols, hdrs] = cbu_dicom_convert(dicom_directory, output_directory, filt, output_format)

if nargin < 1
  dicom_directory = [];
end
if isempty(dicom_directory)
  dicom_directory = spm_select(1, 'dir', 'DICOM directory');
end
if nargin < 2
  output_directory = [];
end
if isempty(output_directory)
  output_directory = pwd;
end
output_directory = spm_select('CPath', output_directory);
if nargin < 3
  filt = [];
end
if isempty(filt)
  filt = '^.*\.dcm$';
end
if nargin < 4
  output_format = [];
end
if isempty(output_format)
  output_format = 'img';
end
cwd = pwd;
output_filt = ['^.*\.' output_format '$'];
files= cellstr(spm_select('List', dicom_directory, filt));
if isempty(files)
  return
end
% Adding full path name to header filenames
for fn = 1:length(files)
  files{fn} = fullfile(dicom_directory, files{fn});
end
hdrs = spm_dicom_headers(char(files));
% SPM dicom converter does not return what files it has written
before_convert = dir(output_directory);
if ~exist(output_directory, 'dir')
  [p f e] = fileparts(output_directory);
  if isempty(p)
    p = pwd;
  end
  mkdir(output_directory)
end
cd(output_directory);
spm_dicom_convert(hdrs, 'all', 'flat', output_format);
cd(cwd);
after_convert = dir(output_directory);
new_files = sf_newer_files(after_convert, before_convert, output_filt);
filenames = char(strcat(output_directory, filesep, new_files));
% Write DTI information if it appears to be present
if sf_is_siemens_dwi(hdrs{1})
  [bvals bvecs] = cbu_dti_params(hdrs);
  prefix = deblank(filenames(1,:));
  [p prefix e] = fileparts(prefix);
  cbu_write_fdt(bvals, bvecs, output_directory, [prefix '.']);
end
return

function new_files = sf_newer_files(new_dirstruct, old_dirstruct, filt)
% return files matching filt in written_path, newer than prior_dir list
new_dirs_filtered = sf_filter_dir(new_dirstruct, filt);
old_dirs_filtered = sf_filter_dir(old_dirstruct, filt);
new_files = {};
if isempty(new_dirs_filtered)
  return
end
new_names = {new_dirs_filtered.name};
if isempty(old_dirs_filtered)
  old_names = {};
else
  old_names = {old_dirs_filtered.name};
end
for f = 1:length(new_dirs_filtered)
  old_element = find(strcmp(new_names{f}, old_names));
  if ~isempty(old_element)
    if datenum(new_dirs_filtered(f).date) <= ...
	  datenum(old_dirs_filtered(old_element).date)
      continue
    end
  end
  new_files{end+1} = new_names{f};
end
return

function filtered_dir = sf_filter_dir(dir_struct, filt)
names = {dir_struct.name};
filtered_dir = [];
for nn = 1:length(names)
  if regexp(names{nn}, filt)
    filtered_dir = [filtered_dir dir_struct(nn)];
  end
end
return

function tf = sf_is_siemens_dwi(hdr)
% returns 1 if this appears to be Siemens DWI sequence
try
  tf = strcmp(hdr.CSAImageHeaderInfo(78).name, 'B_matrix');
catch
  tf = 0;
end
return

