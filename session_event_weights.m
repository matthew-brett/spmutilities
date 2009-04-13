function sew = session_event_weights(sess_event_numbers)
% Calculate equivalent weights given number of events per session
  
n_sessions = size(sess_event_numbers, 1);
n_events = sum(sess_event_numbers, 1);
wt = repmat(n_events, n_sessions, 1);
rows_with_events = logical(n_events ~=0);
sew = zeros(size(sess_event_numbers));
sew(:, rows_with_events) = sess_event_numbers(:, rows_with_events) ./ ...
      wt(:, rows_with_events) * n_sessions;
return
