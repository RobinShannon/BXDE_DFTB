#!/bin/bash
# request resources:
#PBS -l nodes=1:ppn=1
#PBS -l walltime=200:00:00
# on compute node, change directory to 'submission directory':
cd $PBS_O_WORKDIR
cp tmp dftb_in.hsd 
cp dftb_in.hsd tmp
# run your program, timing it for good measure:
i="0"
rm qsub.sh.*
#Path for output
pa=$(echo "/panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/dftb+_1.2.2_autotest/md/IsoRes/1500Res/Paths/")
export OMP_NUM_THREADS=1
while [ $i -lt 100 ]
do

rm Frag*
rm newMin.xyz
EneF1=0
EneF2=0
EneF3=0
EneTS=0

../../../prg_dftb/_obj_x86_64-linux-ifort/dftb+
i=$[$i+1]
cp geom.out.xyz dyn.xyz

Steps=$(awk '/Step:/{ print $2};' transition.xyz)

#Split xyz output into disconnected frgments
/cm/shared/apps/OpenBabel-2.3.2/bin/babel -ixyz newMin.xyz --separate -m -oxyz Frag.xyz

#Convert newMin to gen format
 ../dptools-0.1p3/build/scripts-2.7/xyz2gen newMin.xyz

 
#Relax hygrogens at product geom

if [ -f "Frag2.xyz" ];
then
   cat ../IsoRes/Scripts/GeomTop.txt newMin.gen ../IsoRes/Scripts/GeomAllBottom.txt > geoTemp
else
   cat ../IsoRes/Scripts/GeomTop.txt newMin.gen ../IsoRes/Scripts/GeomAllBottom.txt > geoTemp
fi
cp geoTemp dftb_in.hsd  
../../../prg_dftb/_obj_x86_64-linux-ifort/dftb+ dftb_geom_in.hsd
EneNon=$(awk '/Total energy:/{ print $3};' detailed.out)
cp geom.out.xyz partialOpt.xyz

#Convert Fragments to Gen format and minimise
#One
 ../dptools-0.1p3/build/scripts-2.7/xyz2gen Frag1.xyz
cat ../IsoRes/Scripts/GeomTop.txt Frag1.gen ../IsoRes/Scripts/GeomAllBottom.txt > geoTemp
cp geoTemp dftb_in.hsd
../../../prg_dftb/_obj_x86_64-linux-ifort/dftb+ dftb_geom_in.hsd
EneF1=$(awk '/Total energy:/{ print $3};' detailed.out)
cp geom.out.xyz Frag1Min.xyz
Number1=$(awk 'NR==1{print ($1*3) - 6}' Frag1Min.xyz)
cp Hessian.inp dftb_in.hsd  
../../../prg_dftb/_obj_x86_64-linux-ifort/dftb+
VibF1=$(/panfs/panasas01/chem/rs15164/modes_1.2.2/modes_1.2.2.x86_64-linux | awk -v var="$Number1" '/Vibrational/ {for(i=1; i<=var; i++) {getline; print}}' | tail -n +7)
/cm/shared/apps/OpenBabel-2.3.2/bin/babel -ixyz Frag1Min.xyz -ocml | tail -n +3 | head -n -1 > xmlGeomF1
xmlName1=$(/cm/shared/apps/OpenBabel-2.3.2/bin/babel -ixyz Frag1Min.xyz -ocan | awk '{print $1}' | tr --delete '\\'  | tr --delete '\/' | tr --delete '@' | tr '\(' 'b' | tr '\)' 'b')

#two
if [ -f "Frag2.xyz" ];
then
 ../dptools-0.1p3/build/scripts-2.7/xyz2gen Frag2.xyz
cat ../IsoRes/Scripts/GeomTop.txt Frag2.gen ../IsoRes/Scripts/GeomAllBottom.txt > geoTemp
cp geoTemp dftb_in.hsd
../../../prg_dftb/_obj_x86_64-linux-ifort/dftb+ dftb_geom_in.hsd
EneF2=$(awk '/Total energy:/{ print $3};' detailed.out)
cp geom.out.xyz Frag2Min.xyz
Number2=$(awk 'NR==1{print ($1*3) - 6}' Frag2Min.xyz)
cp Hessian.inp dftb_in.hsd  
../../../prg_dftb/_obj_x86_64-linux-ifort/dftb+
if [Number2 != 1]
then
VibF2=$(/panfs/panasas01/chem/rs15164/modes_1.2.2/modes_1.2.2.x86_64-linux | awk -v var="$Number2" '/Vibrational/ {for(i=1; i<=var; i++) {getline; print}}' | tail -n +7)
fi
/cm/shared/apps/OpenBabel-2.3.2/bin/babel -ixyz Frag2Min.xyz -ocml | tail -n +3 | head -n -1 > xmlGeomF2
xmlName2=$(/cm/shared/apps/OpenBabel-2.3.2/bin/babel -ixyz Frag2Min.xyz -ocan | awk '{print $1}' | tr --delete '\\'  | tr --delete '\/' | tr --delete '@' | tr '\(' 'b' | tr '\)' 'b')
fi

#three
if [ -f "Frag3.xyz" ];
then
 ../dptools-0.1p3/build/scripts-2.7/xyz2gen Frag3.xyz
cat ../IsoRes/Scripts/GeomTop.txt Frag3.gen ../IsoRes/Scripts/GeomAllBottom.txt > geoTemp
cp geoTemp dftb_in.hsd
../../../prg_dftb/_obj_x86_64-linux-ifort/dftb+ dftb_geom_in.hsd
EneF3=$(awk '/Total energy:/{ print $3};' detailed.out)
cp geom.out.xyz Frag3Min.xyz
cp Hessian.inp dftb_in.hsd  
../../../prg_dftb/_obj_x86_64-linux-ifort/dftb+
VibF3=$(/panfs/panasas01/chem/rs15164/modes_1.2.2/modes_1.2.2.x86_64-linux | awk '/Vibrational/ {for(i=1; i<=51; i++) {getline; print}}' | tail -n +7)
xmlGeomF3=$(/cm/shared/apps/OpenBabel-2.3.2/bin/babel -ixyz Frag3Min.xyz -ocml | tail -n +3 | head -n -1)
xmlName3=$(/cm/shared/apps/OpenBabel-2.3.2/bin/babel -ixyz Frag3Min.xyz -ocan | awk '{print $1}' | tr --delete '\\'  | tr --delete '\/' | tr --delete '@' | tr '\(' 'b' | tr '\)' 'b')
fi
  

#Write Smiles
smi1=$(/cm/shared/apps/OpenBabel-2.3.2/bin/babel -ixyz Frag1Min.xyz -ocan | awk '{print $1}' | tr --delete '\\'  | tr --delete '\/' | tr --delete '@' | tr '\(' 'b' | tr '\)' 'b')
if [ -f "Frag2.xyz" ];
then
smi2=$(/cm/shared/apps/OpenBabel-2.3.2/bin/babel -ixyz Frag2Min.xyz -ocan | awk '{print $1}' | tr --delete '\\'  | tr --delete '\/' | tr --delete '@' | tr '\(' 'b' | tr '\)' 'b')
diss=$(echo "diss")
smile2=$diss$smi2
fi
if [ -f "Frag3.xyz" ];
then
smi3=$(/cm/shared/apps/OpenBabel-2.3.2/bin/babel -ixyz Frag3Min.xyz -ocan | awk '{print $1}' | tr --delete '\\'  | tr --delete '\/' | tr --delete '@' | tr '\(' 'b' | tr '\)' 'b')
smile3=$diss$smi3
fi

#Get new directory name
if [ -f "Frag2Min.xyz" ];
then
   if [ -f "Frag3Min.xyz" ];
   then
      title="$pa$smi1$smile2$smile3"
      sm="$smi1$smile2$smile3" 
   else
      title="$pa$smi1$smile2"
      sm="$smi1$smile2"
   fi 
else
   title="$pa$smi1"
   sm="$smi1"
fi

if [ -d "$title" ]
   then
   echo -e "$Steps\t$EneNon\t$EneF1\t$EneF2\t$EneF3" >> $title/Summary.txt
else
   mkdir "$title"
   cat partialOpt.xyz >> $title/minGuess
   cat Frag1Min.xyz >> $title/Frag1.xyz
   cat Frag2Min.xyz >> $title/Frag2.xyz
   cat Frag3Min.xyz >> $title/Frag3.xyz
   cat dyn.xyz >> $title/dyn.xyz	
   cat ../IsoRes/Scripts/MoleT.xml xmlGeomF1 ../IsoRes/Scripts/MoleB.xml > $title/Frag1.xml
   perl -i -pe 's/VIB/'"$VibF1"'/g' $title/Frag1.xml
   sed -i -e "s/NAME/$xmlName1/g" $title/Frag1.xml
   EneF1=$(echo $EneF1" 2625.5" | awk '{printf "%4.3f\n",$1*$2}')
   sed -i -e "s/ENE/$EneF1/g" $title/Frag1.xml
   cat $title/Frag1.xml >> $pa/Molecules.xml
   if [ -f "Frag2.xyz" ];
   then
   cat ../IsoRes/Scripts/MoleT.xml xmlGeomF2 ../IsoRes/Scripts/MoleB.xml > $title/Frag2.xml
   sed -i -e 's/NAME/'$xmlName2'/g' $title/Frag2.xml
   perl -i -pe 's/VIB/'"$VibF2"'/g' $title/Frag2.xml
   EneF2=$(echo $EneF2"  2625.5" | awk '{printf "%4.3f\n",$1*$2}')
   EneF3=$(echo $EneF3"  2625.5" | awk '{printf "%4.3f\n",$1*$2}')
   EneF2=$(echo $EneF2 + $EneF3 | bc)
   sed -i -e 's/ENE/'$EneF2'/g' $title/Frag2.xml
   cat $title/Frag2.xml >> $pa/Molecules.xml
   cp /panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/dftb+_1.2.2_autotest/md/IsoRes/Scripts/Reac2.xml $title/Reac.xml
   sed -i -e 's/PROD1/'$xmlName1'/g' $title/Reac.xml
   sed -i -e 's/PROD2/'$xmlName2'/g' $title/Reac.xml
   sed -i -e 's/ii/'$i'/g' $title/Reac.xml
   else
   cp /panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/dftb+_1.2.2_autotest/md/IsoRes/Scripts/Reac1.xml $title/Reac.xml
   sed -i -e 's/PROD1/'$xmlName1'/g' $title/Reac.xml
   sed -i -e 's/ii/'$i'/g' $title/Reac.xml
   fi
   cat transition.xyz >> $title/transGuess.xyz
   echo -e "$Steps\t$EneNon\t$EneF1\t$EneF2\t$EneF3" >> $title/Summary.txt
   total=$(echo $EneF1 + $EneF2 | bc)
   # Transition state calc
   sed -i '2s/.*/12.0 12.0 1.0 12.0 1.0 1.0 12.0 1.0 1.0 1.0 12.0 1.0 1.0 16.0 1.0 16.0 16.0/' partialOpt.xyz
   cp partialOpt.xyz /panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/dftb+_1.2.2_autotest/md/IsoRes/Transition
   cd /panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/dftb+_1.2.2_autotest/md/IsoRes/Transition
   /panfs/panasas01/chem/rs15164/qct/source/x_Lynx -f path.inp 
   wait    
   EneTS=$(awk '/Transition/{nr[NR+2]}; NR in nr' path.log )
   EneProd=$(awk '/Transition/{nr[NR+606]}; NR in nr' path.log )
   if [ "$EneTS" != "$EneProd" ];
   then
      awk '/Transition/ {for(i=1; i<=19; i++) {getline; print}}' path.log | cat >> /panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/dftb+_1.2.2_autotest/md/1500Res/Paths/$sm/tsMin.xyz
      awk '/Transition/ {for(i=1; i<=19; i++) {getline; print}}' path.log | cat > /panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/dftb+_1.2.2_autotest/md/Isoprene1000bxd/tsMin.xyz
   fi
   cp path.log /panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/dftb+_1.2.2_autotest/md/1500Res/Paths/$sm
   cd /panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/dftb+_1.2.2_autotest/md/Isoprene1000bxd
   echo -e "$Steps\t$sm\t$total\t$EneTS" >> $pa/Energies.txt
   if [ -f "$title/tsMin.xyz" ];
   then
       ../dptools-0.1p3/build/scripts-2.7/xyz2gen tsMin.xyz
       cat ../IsoRes/Scripts/GeomTop.txt tsMin.gen ../IsoRes/Scripts/GeomAllBottom.txt > geoTemp
       cp tsMin.gen geom.out.gen
       cp Hessian.inp dftb_in.hsd  
       ../../../prg_dftb/_obj_x86_64-linux-ifort/dftb+
       VibTS=$(/panfs/panasas01/chem/rs15164/modes_1.2.2/modes_1.2.2.x86_64-linux | awk '/Vibrational/ {for(i=1; i<=51; i++) {getline; print}}' | tail -n +8)
       VibIM=$(/panfs/panasas01/chem/rs15164/modes_1.2.2/modes_1.2.2.x86_64-linux | awk '/Vibrational/ {for(i=1; i<=1; i++) {getline; print}}')
       rm xmlGeomTS 
       /cm/shared/apps/OpenBabel-2.3.2/bin/babel -ixyz tsMin.xyz -ocml | tail -n +3 | head -n -1 > xmlGeomTS
       cat /panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/dftb+_1.2.2_autotest/md/IsoRes/Scripts/TST.xml xmlGeomTS ../IsoRes/Scripts/TSB.xml > $title/TS.xml
       sed -i -e "s/NAME/$sm/g" $title/TS.xml
       perl -i -pe 's/VIB/'"$VibTS"'/g' $title/TS.xml
       sed -i -e "s/IM/$VibIM/g" $title/TS.xml
       sed -i -e "s/ENE/$EneTS/g" $title/TS.xml
       cat $title/TS.xml >> $pa/Molecules.xml
       sed -i -e "s/Trans/$sm/g" $title/Reac.xml
       if [ -f "$title/Frag2.xyz" ];
       then
           sed -i -e '11,16d;' $title/Reac.xml
       fi
   else
       if [ -f "$title/Frag2.xyz" ];
       then
           sed -i -e '17,20d;' $title/Reac.xml
       fi             
   fi
   cat $title/Reac.xml >> $pa/Reactions.xml 
fi
cd /panfs/panasas01/chem/rs15164/DFTB+/dftb+_1.2.2_src/dftb+_1.2.2_autotest/md/Isoprene1000bxd
cp tmp dftb_in.hsd
done
