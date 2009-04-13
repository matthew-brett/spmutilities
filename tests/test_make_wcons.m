% Script to show weighted contrast stuff
n_conds = 24; % 24 conditions
n_extras = 6;
fname = data_path('sesscon_txt_example.txt');
sess_cons = load_sesscon_txt(fname, n_conds);
ons_fname = data_path('sesscon_n_ons');
load(ons_fname);
sess_ev_weights = session_event_weights(n_ons);
weighted_contrasts = weight_contrasts(sess_cons, ...
				      sess_ev_weights, ...
				      n_extras);
% Unequal number of extra regressors per session
n_extras = [6 6 9 10 6];
wc2 = weight_contrasts(sess_cons, ...
		       sess_ev_weights, ...
		       n_extras);
