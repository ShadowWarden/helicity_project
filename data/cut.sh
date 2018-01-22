#!/bin/bash

gtselect infile=photon_files.txt \
	evclass=128 evtype=3\
	outfile=gtselected_P8_50-60GeV.fits \
	ra=INDEF dec=INDEF rad=INDEF \
	tmin=INDEF tmax=INDEF \
	emin=50e3 emax=60e3 zmax=180
