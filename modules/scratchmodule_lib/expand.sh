#!/bin/bash 

filename_rules=$1
filename_templ=$2
fileout=${filename_rules/rules.sh/inc}
rm -f $fileout
source ${filename_rules} | while read line; do
    eval "$line"
    cat ${filename_templ} >> $fileout
    for key in `echo "$line" | awk -F= '{for (i=1; i<NF; i++) {print $i}}' | awk '{print $NF}'`;
    do
	sed -i 's/${'${key}'}/'"${!key}"'/g' $fileout
    done
done 
echo "SUCCESS: The generated files is called" $fileout "!"

