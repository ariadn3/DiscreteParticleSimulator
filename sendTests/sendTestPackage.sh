scp -r -q -o LogLevel=QUIET $1 $2:~/$1
ssh -q $2 "cd ~/$1; nohup ./run.sh > /dev/null &"
