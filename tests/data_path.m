function dpath = data_path(fname)
% Return path to test data, optionally with filename
if nargin < 1
  fname = [];
end
my_path = fileparts(mfilename('fullpath'));
dpath = fullfile(my_path, 'data');
if ~isempty(fname)
  dpath = fullfile(dpath, fname);
end
return

