#!/bin/sh
function incar_xc(){
if [ "$1" = PE ] ;then 
cat > pbe <<qwer
#PBE
qwer
cat pbe >>INCAR ; rm pbe
elif [ "$1" = PS ] ; then
cat > pbe <<qwer
#PBEsol
GGA=PS
qwer
cat pbe >>INCAR ; rm pbe
elif [ "$1" = revTPSS ] ; then
    echo "metaGGA is CHOSEN "$1
    if [[ `grep " kinetic energy-density" POTCAR` ]]; then
    cat > pbe <<qwer
# revTPSS
METAGGA=RTPSS  
qwer
cat pbe >>INCAR ; rm pbe
    else
        echo "WARNING !!!! THE FUNCTIONAL is NOT CHOSEN CORRECTLY"
        echo "metaGGA is CHOSEN but there is no kinetic information"
        echo
        echo "PLEASE CHECK THE INPUT incar_xc "
        echo "PLEASE CHECK THE POTCAR, there is no kinetic information"
        #exit 0
    fi
 
elif [ "$1" = lda ]; then
echo LDA is chosen
    if [[ `grep " PAW " POTCAR` ]]; then
        echo "POTCAR for LDA is chosen CORRECTLY"
    else
        echo "POTCAR for LDA is not chosen CORRECTLY, please check again"
    fi
fi
 
}
 
function incar_vdw(){
# $1 = specific number of vdW correction
# $2 = path of vasp_unitiliy including vdw_kernel.bindat 
if [ "$1" = 10 ] ;then 
cat > vdW <<qwer
#DFT-D2
IVDW=10
qwer
vdw=D2
elif [ "$1" = 11 ] ; then
cat > vdW <<qwer
#DFT-D3_GRIMME
IVDW=11
qwer
vdw=D3
elif [[ "$1" = 12 ]]; then
cat > vdW <<qwer
#DFT-D3_BeckeJonson
IVDW=12
qwer
vdw=D3BeckeJonson
elif [[ "$1" = 20 ]]; then
cat > vdW <<qwer
#vdW_TS
IVDW=20
qwer
vdw=TS
elif [[ "$1" = 21 ]]; then
cat > vdW <<qwer
# Van der Waals of TS-SCS
IVDW=2
LVDWSCS = .TRUE.
qwer
vdw=TS-SCS
elif [[ "$1" = 30 ]]; then
cat > vdW <<qwer
#optB86b-vdW (2011)
GGA=MK
PARAM1 = 0.1234
PARAM2 = 1.0000
LUSE_VDW = .TRUE.
AGGAC = 0.0000  ## 
qwer
vdw=optB86b
cp $2/vdw_kernel.bindat . 
elif [[ "$1" = 31 ]]; then
cat > vdW <<qwer
#optB88-vdW 
GGA=BO
PARAM1 = 0.1833333333
PARAM2 = 0.2200000000
LUSE_VDW = .TRUE.
AGGAC = 0.0000
qwer
vdw=optB88
cp $2/vdw_kernel.bindat . 
elif [[ "$1" = 32 ]]; then
cat > vdW <<qwer
#optPBE-vdW (DFT-DF vdW)
GGA=OR
LUSE_VDW = .TRUE.
AGGAC = 0.0000
qwer
vdw=optPBE
cp $2/vdw_kernel.bindat . 
elif [[ "$1" = 33 ]]; then
cat > vdW <<qwer
#DFT-DF vdW (revPBE)   
GGA=RE
LUSE_VDW = .TRUE.
AGGAC = 0.0000
qwer
vdw=vdW_DF1_revPBE
cp $2/vdw_kernel.bindat . 
elif [[ "$1" = 34 ]]; then
cat > vdW <<qwer
#vdW-DF2
GGA=ML
LUSE_VDW = .TRUE.
Zab_vdW = -1.8867
AGGAC = 0.0000
qwer
vdw=vdW_DF2
cp $2/vdw_kernel.bindat . 
else
### ADDed 
elif [[ "$1" = 300 ]]; then
cat > vdW <<qwer
#optB86b-vdW (2011)
GGA=MK
# PARAM1 = 0.1234
# PARAM2 = 1.0000
# LUSE_VDW = .TRUE.
# AGGAC = 0.0000  ## 
qwer
vdw=optB86b
elif [[ "$1" = 310 ]]; then
cat > vdW <<qwer
#optB88-vdW 
GGA=BO
# PARAM1 = 0.1833333333
# PARAM2 = 0.2200000000
# LUSE_VDW = .TRUE.
# AGGAC = 0.0000
qwer
vdw=optB88
elif [[ "$1" = 320 ]]; then
cat > vdW <<qwer
#optPBE-vdW (DFT-DF vdW)
GGA=OR
# LUSE_VDW = .TRUE.
# AGGAC = 0.0000
qwer
vdw=optPBE
elif [[ "$1" = 330 ]]; then
cat > vdW <<qwer
#DFT-DF vdW (revPBE)   
GGA=RE
LUSE_VDW = .TRUE.
AGGAC = 0.0000
qwer
vdw=vdW_DF1_revPBE
elif [[ "$1" = 340 ]]; then
cat > vdW <<qwer
#vdW-DF2
GGA=ML
LUSE_VDW = .TRUE.
Zab_vdW = -1.8867
AGGAC = 0.0000
qwer
vdw=vdW_DF2
else
 
vdw=pure
cat > vdW <<qwer
#pure
qwer
fi
cat vdW #>>INCAR ; rm vdW
 
}
 
 
function incar_atom(){
cat > ./INCAR <<qwer
PREC=accurate
ENCUT=500
NELMIN=7     
EDIFF=1E-5
NSW=0       ##
IBRION=-1      ##relaxation allowed
ADDGRID=T
LREAL=F
LASPH=T
LORBIT=10
 
##Determining the groundstate energy of atoms:
ISPIN=2
ISYM=0 #no symmetry
ISMEAR = 0
SIGMA=0.002
AMIX=0.2
BMIX=0.0001
NELM=100
ICHARG=1
 
qwer
}
 
 
 
function incar_vcopt(){
cat > ./INCAR <<qwer
PREC=accurate
ENCUT=500
NELMIN=5
EDIFF=1E-5
NSW=150
IBRION=2
POTIM=0.3
ISIF=3
ISMEAR=1
SIGMA=0.1
#ADDGRID=T
LWAVE=F
LCHARG=F
LREAL=F
LASPH=T
LORBIT=10
 
 
qwer
}
function incar_scf(){
cat > ./INCAR <<qwer
PREC=accurate
ENCUT=500
NELMIN=10
EDIFF=1E-5
NSW=0
IBRION=-1
ISMEAR=1
SIGMA=0.1
ADDGRID=T
LWAVE=F
LCHARG=F
LREAL=F
LASPH=T
LORBIT=10
 
qwer
}
 
function incar_scf_tet(){
cat > ./INCAR <<qwer
PREC=accurate
ENCUT=500
NELMIN=7
EDIFF=1E-5
NSW=0
IBRION=-1
ISMEAR=-5  #-4-tet -1-fermi 0-gaus
SIGMA=0.01 #broadening in eV
#ADDGRID=T
LWAVE=F   #write WAVECAR
LCHARG=F  #write CHGCAR
LORBIT=10
LASPH=T
LREAL=F
 
qwer
}
 
function incar_dos(){
cat > ./INCAR <<qwer
PREC=accurate
ENCUT=500
NELMIN=10
EDIFF=1E-5
NSW=0
IBRION=-1
ISMEAR=-5
SIGMA=0.01
#ADDGRID=T
LWAVE=F
LCHARG=F
LREAL=F
LASPH=T
LORBIT=10
 
EMIN=-20.00
EMAX=20.00
NEDOS=501
qwer
}
 
 
function incar_ev(){
cat > ./INCAR <<qwer
PREC=accurate
ENCUT=500
NELMIN=10
EDIFF=1E-5
NSW=150
IBRION=2
POTIM=0.1   # MTGstandard=0.5
ISIF=4
ISMEAR=1
SIGMA=0.1
ADDGRID=T
LWAVE=F
LCHARG=F
LREAL=F
LASPH=T
LORBIT=10
 
qwer
}
 
 
function incar_stiff(){
cat > ./INCAR <<qwer
PREC=accurate
ENCUT=500
NELMIN=5
EDIFF=1E-5
NSW=1
IBRION=6
POTIM=0.01
ISIF=3
ISMEAR=-5
SIGMA=0.001
#ADDGRID=T
LWAVE=F
LCHARG=F
LREAL=F
LASPH=T
LORBIT=10
NFREE=4
 
qwer
}
 
 
function incar_geoopt(){
cat > ./INCAR <<qwer
PREC=accurate
ENCUT=500
NELMIN=5
EDIFF=1E-5
NSW=150
IBRION=2
POTIM=0.1
ISMEAR=1
SIGMA=0.1
#ADDGRID=T
LWAVE=F
LCHARG=F
LREAL=F
LASPH=T
LORBIT=10
 
qwer
}
