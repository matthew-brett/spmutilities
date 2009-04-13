fname = data_path('sesscon_txt_example.txt');
n_conds = 20;
sess_cons = load_sesscon_txt(fname, n_conds);
assert_equal(size(sess_cons), [1 25]);
eg_con = sess_cons(3);
assert_strings_equal(eg_con.name, 'IDMiss_NC');
assert_equal(numel(eg_con.contrast), n_conds);




