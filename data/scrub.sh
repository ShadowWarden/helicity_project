#!/bin/bash

gtmktime spacecraft_files.txt \
	filter="DATA_QUAL>0 && LAT CONFIG==1" \
	roicut=no \
	evfile=gtselected_P8_10-20GeV.fits \
	outfile=FermiData_PASS8_10-20GeV.fits
