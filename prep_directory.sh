#!/bin/bash

# Run this script to prepare the orthoCapture directory for running orthoCapture

NTLOC=$1
./ntdb_location.sh ${NTLOC}

./uncompress_blast.sh
	# uncompresses blast files
