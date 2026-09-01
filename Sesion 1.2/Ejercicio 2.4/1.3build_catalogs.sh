#!/bin/bash

src=$1
files=$2

# 
# Use a loop to create a list of files to supply to cstacks.
# 
samp="" 
for file in $files 
do 
    samp+="-s $src/stacks/$file "; 
done 

# 
# Build the catalog; the "&>>" will capture all output and append it to the Log file.
# 
cstacks -g -p 36 -b 1 -n 1 -o $src/stacks $samp &>> $src/stacks/Log 
