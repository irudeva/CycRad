:
f90 -o trk_all.exe trk_all.f
HS=NH
indir="/projects/drakkar/cyclone/Track/Original/${HS}old/"
outdir="/projects/drakkar/cyclone/Track/Directold/"
:
firstyear=2005
lastyear=2007
year=${firstyear}
while [ $year -le $lastyear ]
do
echo $year
#if ls ${indir} | grep int${year}-up.trk.gz
if ls ${indir} | grep $HS.${year}.trk.gz
then
#gunzip ${indir}int${year}-up.trk.gz
gunzip ${indir}$HS.${year}.trk.gz
fi

trk_all.exe<<mark
${indir}$HS.${year}.trk
${outdir}num${year}.$HS.d
${outdir}${year}.$HS.drk
mark
:
gzip ${indir}$HS.${year}.trk
: 
year=` expr $year + 1 `
done
rm trk_all.exe
