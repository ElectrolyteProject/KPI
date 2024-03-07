# A script to calculate RESP or RESP2 charges by invoking Multiwfn
# Written by Tian Lu (sobereva@sina.com), 2020-Jun-30
# Modified by Nan Yao, 2021-Dec-17
# Example of calculation of RESP charges: ./calcRESP.sh H2CO_gas.fchk
# Example of calculation of RESP2 charges: ./calcRESP.sh H2CO_gas.fchk H2CO_water.fchk

#!/bin/bash
echo Calculating RESP charge for $1 ...
dirfchk=$2/$1
Multiwfn $dirfchk << EOF > /dev/null
7
18
1
y
0
0
q
EOF
chgfile=${1%.*}.chg
dirchgfile=$2/$chgfile

if [ $3 ];then
#echo
echo Calculating RESP charge for $3 ...
dirfchk2=$4/$3
Multiwfn $dirfchk2 << EOF > /dev/null
7
18
1
y
0
0
q
EOF
chgfile2=${3%.*}.chg
dirchgfile2=$4/$chgfile2

delta=0.5
if [ $5 ];then
	delta=$5
fi
#echo
echo delta parameter is $delta

dirresp2=$4/RESP2.chg

#Calculate RESP2 charge
paste $chgfile $chgfile2 |awk '{printf $1 " " $2 " " $3 " " $4 " " (1-d)*$5+d*$10 "\n"}' d=$delta > $dirresp2
echo RESP2 charges has been written to RESP2.chg

#mv $chgfile $dirchgfile
#mv $chgfile2 $dirchgfile2
#echo RESP charges for $1 and $3 have been moved to corresponding directory
fi
