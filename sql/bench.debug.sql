select i_lev_exc,i_lev_test_real,avg(taxid_lev_test==taxid_lev_pred) as avg_acu from bench_samp group by i_lev_exc,i_lev_test_real;
select i_lev_exc,i_lev_pred,avg(taxid_lev_test==taxid_lev_pred) as avg_acu from bench_samp group by i_lev_exc,            i_lev_per;

