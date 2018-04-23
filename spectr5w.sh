#!/bin/bash
export SPECTR_NSECT=1
export SPECTR_BGSCALE=2.24
export SPECTR_AUXCUT=""
export SPECTR_OUTDIR=period_all
export SPECTR_PAIRDIR="/mnt/root1/danss_pair6/"
for ((i=1; $i<235; i=$i+1 )) ; do
	./spectr5w $i
done

