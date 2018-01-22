#!/bin/bash

A=151
B=200
download="/home/omkar/programs/helicity/data/photon/"

for((i=$A; i<=$B; ++i))
do 
	a="wget -P $download https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/weekly/photon/lat_photon_weekly_w"
	b=$i"_p302_v001.fits"
	$a$b
done
