# ~/nci_gnu_parallel 'Rscript -e "print(\"hi\")"'

if [ "$1" = "" ]; then echo -e "USAGE:\nRun this command first: for i in \$(compgen -v); do export \$i; done\nARGUMENTS:\n. Wildcards must always have a quote around them: e.g. \"xxx\"\*\".txt\". Escapes required for: \$, \", \\, \*.\n\$1: total no. processes running at once\n\$2: no. processes per node\n\$3: no. cores exclusively assigned to each job\n\$4: no. chunks\n"; exit 2; fi

echo \#\!/bin/bash > $TMPDIR"/nci_multinode_envfile.txt"

# printenv | sed 's|^\([^=]*\)=\(.*\)|\1="\2"|g' >> $TMPDIR"/nci_multinode_envfile.txt"

printenv | grep ^PATH= | sed 's|^\([^=]*\)=\(.*\)|export \1="\2"|g' >> $TMPDIR"/nci_multinode_envfile.txt"

printenv | grep ^LIBRARY_PATH= | sed 's|^\([^=]*\)=\(.*\)|export \1="\2"|g' >> $TMPDIR"/nci_multinode_envfile.txt"

printenv | grep ^LD_LIBRARY_PATH= | sed 's|^\([^=]*\)=\(.*\)|export \1="\2"|g' >> $TMPDIR"/nci_multinode_envfile.txt"

printenv | grep ^TMPDIR= | sed 's|^\([^=]*\)=\(.*\)|export \1="\2"|g' >> $TMPDIR"/nci_multinode_envfile.txt"

# either force all commands to be on one line
# OR
# we need to replace all backslashes with double backslashes 
read -e -p "GNU Parallel command here:
" cmd_gnu_parallel

echo $cmd_gnu_parallel > $TMPDIR"tmp_command0.txt"

echo "====> Input command:"
cat $TMPDIR"tmp_command0.txt"
echo -e "\n"

cat $TMPDIR"tmp_command0.txt" | sed "s|^parallel|parallel --dry-run|g" > $TMPDIR"tmp_command1.txt"

echo "====> Dry run command:"
cat $TMPDIR"tmp_command1.txt"
echo -e "\n"

# user verify that commands parsed correctly
while true; do
    read -e -p "Confirm that commands have passed correctly (Y/N): " stdin_user_verify
    if [ "$stdin_user_verify" = "Y" ]; then
        echo -e "\nExecuting...\n"
        break
    elif [ "$stdin_user_verify" = "N" ]; then
        exit
    else
        echo -e "\nType in Y or N\n"
    fi
done

cat $TMPDIR"/nci_multinode_envfile.txt" $TMPDIR"tmp_command1.txt" > $TMPDIR"tmp_command.txt"

bash $TMPDIR"tmp_command.txt" > $TMPDIR"nci_multinode_cmdfile_L2.txt"

tempfileprefix=$(date '+%s%N')
# $(shuf -i 1-100000 -n 1)

split -a 4 -n l/$4 $TMPDIR"nci_multinode_cmdfile_L2.txt" $TMPDIR"/"$tempfileprefix

for i in $(ls $TMPDIR"/"$tempfileprefix*); do cat $TMPDIR"/nci_multinode_envfile.txt" $i > $i"_withenv"; done

ls $TMPDIR"/"$tempfileprefix*_withenv | awk '{print "bash " $0}' > $TMPDIR"nci_multinode_cmdfile_L1.txt"

mpirun --np $1 --map-by ppr:$2:node:PE=$3 --oversubscribe --bind-to none nci-parallel --input-file $TMPDIR"nci_multinode_cmdfile_L1.txt"

