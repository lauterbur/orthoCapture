#!/bin/bash

j=$1 # directory path
i="'${j}'" # put quotes around directory path
echo '        ntloc = '${i} | sed -i.bak -e '/# Define ntloc/{r /dev/stdin' -e ';p;N;d;}' orthoCapture.py 
        # replaces the line after "# Define ntloc" with ntloc = <newpath>
        # also backs up previous version of orthoCapture to orthoCapture.py.bak

# now orthoCapture knows where your nt database is