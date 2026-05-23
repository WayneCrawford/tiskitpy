#! /bin/csh -f
# L_CHEAPO 2000 FIR filter correction 
# The correction is performed only for the first filter step, 
# which is the only noncausal (?, see Theunisson, RATS cruise
# report 2)
#
# USAGE: corr_lc2000 <input file> <fdig>
# input file: single colummn ASCII format (1 value per line)
# fdig: digitization frequency of the input file in Hz 
# supported values for fdig are: 62.5, 125, 250
#
echo $2 > tmp764995
awk '{printf("%d", $1)}' tmp764995 > tmp956723
set fdig = `cat tmp956723`
if ($2 == 62.5) then
	interpolate -i $1 -o $1.3 -f 2
	fir2caus -c lc2000_fir3_125.prt -i $1.3 -o $1.3c -f $2
	decimate -i $1.3c -o $1.corr -f 0 -d 2
	rm $1.3 $1.3c
else if ($2 == 125) then
	interpolate -i $1 -o $1.3 -f 2
	fir2caus -c lc2000_fir3_250.prt -i $1.3 -o $1.3c -f $2
	decimate -i $1.3c -o $1.corr -f 0 -d 2
	rm $1.3 $1.3c
else
	echo "Unknown sampling rate: $2"
endif
rm tmp764995 tmp956723
