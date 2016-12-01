:
fprog=trk_rad.d2s
f90 -o ${fprog} ${fprog}.f
HS=NH
outdir="/projects/drakkar/cyclone/Track/Rad_real_slp/${HS}old/"
indir="/projects/drakkar/cyclone/Track/Directold/"
:
firstyear=2006
lastyear=2006
year=${firstyear}
while [ $year -le $lastyear ]
do
echo $year
if ls ${indir} | grep rad.${year}.$HS.drk.gz
then
gunzip ${indir}rad.${year}.$HS.drk.gz
fi
if ls ${indir} | grep num${year}.$HS.d.gz
then
gunzip ${indir}num${year}.$HS.d.gz
fi
if ls ${indir} | grep ${year}.$HS.drk.gz
then
gunzip ${indir}${year}.$HS.drk.gz
fi
:
${fprog}<<mark
${indir}num${year}.$HS.d
${indir}${year}.$HS.drk
${indir}rad.${year}.$HS.drk
${outdir}rad.${year}.$HS.trk
mark
:
gzip ${indir}${year}.$HS.drk
gzip ${indir}num${year}.$HS.d
gzip ${indir}rad.${year}.$HS.drk
: 
year=` expr $year + 1 `
done
rm ${fprog}
