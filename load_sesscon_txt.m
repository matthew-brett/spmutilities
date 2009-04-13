function sess_cons = load_sesscon_txt(fname, n_conds)
% Read contrast definition file, return session contrasts
% FORMAT sess_cons = load_sesscon_txt(fname, n_conds)
% 
% Inputs
% fname        - file in format below
% n_conds      - total number of conditions
%
% Outputs
% sess_cons    - contrasts N by n_conds
%
% File should be in format originally used by Eva and Sonia, 
% with lines like:
% my_contrast name     1       8       1       9       1       19
% where the values are tab-separated, the first entry is the contrast
% name (with or without quotes), and the subsequent values are pairs of
% form contrast weight, condition number - e.g 1 8 means apply 1 to
% condition number 8.

fid = fopen(fname);
if fid == -1
  error('Could not read file')
end
con_ind = 1;
while 1
  tline = fgetl(fid);
  if ~ischar(tline), break, end
  tab_inds = strfind(tline, sprintf('\t'));
  if isempty(tab_inds), continue, end
  name = sscanf(tline(1:tab_inds(1)-1), '%s');
  % Remove any quotes from name
  if strcmp(name(1), '''') | strcmp(name(1), '"')
    name = name(2:end-1);
  end
  vals = sscanf(tline(tab_inds(1)+1:end), '%d');
  vals = reshape(vals, 2, length(vals)/2);
  con = zeros(n_conds, 1);
  con(vals(2,:)) = vals(1,:);
  sess_cons(con_ind) = struct(...
      'name',  name,...
      'contrast', con);
  con_ind = con_ind+1;
end
fclose(fid);
return
