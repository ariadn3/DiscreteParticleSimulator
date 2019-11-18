cd $1
rm run.sh
cat ../bashHeader init.sh compile.sh test.sh cleanup.sh > run.sh
