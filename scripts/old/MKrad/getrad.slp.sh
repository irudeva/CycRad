:
ifort -o getrad.slp.exe getrad.slp.f
HS=NH
indir="/projects/drakkar/cyclone/Track/Original/$HS/"
glodir="/projects/drakkar/cyclone/NCEPdata/Slp$HS/Direct/"
outdir="/projects/drakkar/cyclone/Track/Rad_real_slp/$HS/rad."
maskdir="/projects/drakkar/cyclone/181area/"
 
fyr=1948
lyr=1948
yr=${fyr}

while [ $yr -le $lyr ]
do
echo $yr
if [ -e ${glodir}$yr.$HS.dlo.gz ]; then
echo 'gunzip'
gunzip ${glodir}$yr.$HS.dlo.gz; fi
:
if [ -e ${indir}$HS.${yr}.trk.gz ]; then
echo 'gunzip'
gunzip ${indir}$HS.${yr}.trk.gz; fi

getrad.slp.exe <<mark
${indir}$HS.${yr}.trk
${glodir}${yr}.$HS.dlo
${outdir}${yr}.$HS.trk
${maskdir}
mark
gzip ${glodir}$yr.$HS.dlo ${indir}$HS.${yr}.trk
yr=` expr $yr + 1 `
done

rm getrad.slp.exe
