# Runs ssu on the test data and prints out the wanted results.

for f in uwtd1 uwtd2; do
    biom convert --to-hdf5 --table-type="OTU table" -i $f.biom.tsv -o $f.biom
    echo "----- GOT ---------------"
    ssu -m unweighted -i $f.biom -t $f.tree -o /dev/stdout
    echo "----- WANT --------------"
    cat $f.want
done 2> /dev/null

for f in wtd; do
    biom convert --to-hdf5 --table-type="OTU table" -i $f.biom.tsv -o $f.biom
    echo "----- GOT ---------------"
    ssu -m weighted_normalized -i $f.biom -t $f.tree -o /dev/stdout
    echo "----- WANT --------------"
    cat $f.want
done 2> /dev/null

rm *.biom
