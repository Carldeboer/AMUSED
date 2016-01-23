#/bin/bash

AMUSED -s 6 -q test/test1.fa -o test1.res.gz -do
echo "test1: " `zdiff test1.res.gz test/test1.s6.results.gz | wc -l` " differences (should be 0)"

AMUSED -ds -s 6 -q test/test1.fa -o test2.res.gz -do
echo "test2: " `zdiff test2.res.gz test/test2.s6.results.gz | wc -l` " differences (should be 0)"

AMUSED -s 6 -b test/test2.fa.gz  -q test/test1.fa -o test3.res.gz -do
echo "test3: " `zdiff test3.res.gz test/test3.s6.results.gz | wc -l` " differences (should be 0)"

AMUSED -s 6 -b test/test4.fa  -q test/test4_RC.fa -o test4.res.gz -do  -bc -ds
echo "test4: " `zdiff test4.res.gz test/test4.s6.results.gz | wc -l` " differences (should be 0)"

AMUSED-KS -s 6  -q test/test5.fa -o test5.res.gz -ds
echo "test5: " `zdiff test5.res.gz test/test5.s6.results.gz | wc -l` " differences (should be 0)"
