create table test_bulk
(
iid integer,
gi bigint,
taxid integer,
src_db char(1),
project char(4),
seq_len bigint,
acc char(20),
kind char(2),
seq_hdr char(40)
);

copy 50000 offset 0 records into test_bulk from '/home/andrey/work/mgtaxa/test_data/bulk-1.dat';
copy 50000 offset 50000 records into test_bulk from '/home/andrey/work/mgtaxa/test_data/bulk-1.dat';
copy 50000 offset 100000 records into test_bulk from '/home/andrey/work/mgtaxa/test_data/bulk-1.dat';
copy 50000 offset 150000 records into test_bulk from '/home/andrey/work/mgtaxa/test_data/bulk-1.dat';
copy 50000 offset 200000 records into test_bulk from '/home/andrey/work/mgtaxa/test_data/bulk-1.dat';
copy 50000 offset 250000 records into test_bulk from '/home/andrey/work/mgtaxa/test_data/bulk-1.dat';
copy 50000 offset 300000 records into test_bulk from '/home/andrey/work/mgtaxa/test_data/bulk-1.dat';
copy 50000 offset 350000 records into test_bulk from '/home/andrey/work/mgtaxa/test_data/bulk-1.dat';
copy 50000 offset 400000 records into test_bulk from '/home/andrey/work/mgtaxa/test_data/bulk-1.dat';
copy 50000 offset 450000 records into test_bulk from '/home/andrey/work/mgtaxa/test_data/bulk-1.dat';
