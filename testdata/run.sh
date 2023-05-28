# Runs frcfrc on the test data and checks the results.

for f in uwtd1 uwtd2; do
    frcfrc -i $f.dense -t $f.tree -o $f.got
    frcfrc -s -i $f.sparse -t $f.tree -o $f.s.got
done 2> /dev/null

for f in wtd; do
    frcfrc -w -i $f.dense -t $f.tree -o $f.got
    frcfrc -w -s -i $f.sparse -t $f.tree -o $f.s.got
done 2> /dev/null

for f in uwtd1 uwtd2 wtd; do
    diff $f.got $f.want
    diff $f.s.got $f.want
done

rm *.got
