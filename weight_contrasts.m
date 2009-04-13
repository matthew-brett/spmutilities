function weighted_cons = weight_contrasts(sess_cons, sess_wts, n_extras)
% Recalculate contrasts using weights, and omitting absent conditions
% Assumes same (number of) conditions per session
%
% Inputs
% sess_cons     - contrasts to apply for single session assuming all
%                 conditions present. Structure containing field
%                 'contrast', and 'name'.
% sess_wts      - NxM weights to apply, where N is number of sessions
%                 and M is number of conditions (assuming all present 
%                 for every session).  Usually reflects the number of
%                 events in the condition, weighted by the total number
%                 of events in all sessions.  If weights is 0 for a
%                 particular session and condition, that condition is
%                 removed from the contrast vector. 
% n_extras      - scalar, or vector, length N (see above), with number
%                of extra (non-condition) regressors present in this
%                session
%
% Returns
% weighted_cons - contrast, weighted for number of events per condition
%                 in each session (see above)

[n_sessions n_conds] = size(sess_wts);
if numel(n_extras) == 1
  n_extras = ones(1, n_sessions) * n_extras;
end

% NxM weight
wt_absents = sess_wts == 0;
all_wts = sf_to_con(sess_wts, n_extras);
absents = logical(sf_to_con(wt_absents, n_extras));
weighted_cons = sess_cons;
for cno = 1:numel(sess_cons)
  con = sess_cons(cno).contrast';
  con = repmat(con, n_sessions, 1);
  con = sf_to_con(con, n_extras);
  con = con .* all_wts;
  con(absents) = [];
  % Add session effects
  con = [con;zeros(n_sessions,1)];
  if ~any(con)
    % Make a dummy contrast using a block regressor
    con(end) = 1;
    weighted_cons(cno).name = 'Dummy contrast - beware';
  end
  weighted_cons(cno).contrast = con;
end
return

function con = sf_to_con(wts, n_extras)
% Return full contrast vector with zeros for extra regressors
n_sessions = size(wts, 1);
con = [];
for sno = 1:n_sessions
  sess_con = [wts(sno,:), zeros(1, n_extras(sno))];
  con = [con sess_con];
end
con = con';
return
  
