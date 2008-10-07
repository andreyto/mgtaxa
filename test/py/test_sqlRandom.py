"""This was intended to check random sampling w/o replacement with unequal probabilities within SQL.
However, it became obvious that there is no magic in MySQL - select * from tbl order by rand() limit 1
still selects and sorts all records, which is unacceptable for large tables."""

from MGT.Sql import *

db = createDbSql()

db.ddl("create table test_rnd_pop (id int,weight int)",dropList=["table test_rnd_pop"])

db.ddl("""insert into test_rnd_pop
          values
          (1,2),
          (2,1)""")

for i_trial in range(nTrials):
    db.execute("select id from test_rnd_pop order by rand() limit 1")
