:
f90 -o trk_rad.exe trk_rad.f
HS=NH
indir="/projects/drakkar/cyclone/Track/Rad_real_slp/$HS/"
outdir="/projects/drakkar/cyclone/Track/Direct/"
:
firstyear=2005
lastyear=2007
year=${firstyear}
while [ $year -le $lastyear ]
do
echo $year
if ls ${indir} | grep rad.${year}.$HS.trk.gz
then
gunzip ${indir}rad.${year}.$HS.trk.gz
fi

trk_rad.exe<<mark
${indir}rad.${year}.$HS.trk
${outdir}rad.${year}.$HS.drk
mark
:
gzip ${indir}rad.${year}.$HS.trk
: 
year=` expr $year + 1 `
done
rm trk_rad.exe
