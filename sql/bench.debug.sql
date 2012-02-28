select i_lev_exc,i_lev_test_real,avg(taxid_lev_test==taxid_lev_pred) as avg_acu from bench_samp group by i_lev_exc,i_lev_test_real;
select i_lev_exc,i_lev_pred,avg(taxid_lev_test==taxid_lev_pred) as avg_acu from bench_samp group by i_lev_exc,            i_lev_per;
select a.*,b.name from conf_wgt a, taxa_names b where a.taxid_lev_test=b.id and taxid_lev_pred=4890 order by i_lev_exc,i_lev_per,cnt desc;
select a.*,b.name,count(*) as cnt from bench_samp a, taxa_names b where a.taxid_lev_test_bot=b.id and taxid_lev_pred=4890 and taxid_superking <> 2759 group by i_lev_exc,i_lev_per,taxid_lev_test,taxid_lev_test_bot order by i_lev_exc,i_lev_per,cnt desc, taxid_lev_test,taxid_lev_test_bot;

select a.*,b.name,c.sum_cnt as sum_cnt_test,cast(a.cnt as REAL)/c.sum_cnt as ratio_test from conf a, taxa_names b,test_cnt c where a.taxid_lev_test=b.id and a.taxid_lev_test=c.taxid_lev_test and a.i_lev_exc=c.i_lev_exc and a.i_lev_per=c.i_lev_per and taxid_lev_pred=4890 order by a.i_lev_exc,a.i_lev_per,ratio_test desc;

