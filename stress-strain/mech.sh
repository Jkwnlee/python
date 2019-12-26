#!/bin/sh
function FCC_mech(){
#########################################################
##                                                     ##
##  This function is useful for FCC-PRIMITIVE unitcell ##
##  Please use this function with POSCAR and OUTCAR    ##
##  calculated by INCAR including [IBRION=6]           ##
##  Add anisotropic properties of FCC (Elastic, Factor)##
##  Made by Ji-Hwan, Lee                               ##
##  Updated Date: 2015. 2. 22                          ##
#########################################################
# atomE=$( echo $1)
lat_con=0;  c11=0;  c12=0;  c44=0;  s11=0;  s12=0;  s33=0;  v100=0;  v110=0;  v111=0;  A_c=0;  A_z=0;  A_u=0;  
E100=0;  E110=0;  E111=0;  E_E=0;  b_mo_v=0;  b_mo_r=0;  b_mo=0;  y_mo_v=0;  y_mo_r=0;  y_mo=0;  s_mo_v=0;  s_mo_r=0;  s_mo=0;  iso_v=;  
# bulkM=$( echo $2)
egrep -10 ELASTIC OUTCAR | tail -10 >> ./stiffness_report
c11bar=$(head -n 3 ./stiffness_report| tail -n 1 | awk '{print $2}')
c12bar=$(head -n 3 ./stiffness_report| tail -n 1 | awk '{print $3}')
c44bar=$(head -n 6 ./stiffness_report| tail -n 1 | awk '{print $5}')
c11=$( echo "scale=30; ( $c11bar / 10 ) " |bc )
c12=$( echo "scale=30; ( $c12bar / 10 ) " |bc )
c44=$( echo "scale=30; ( $c44bar / 10 ) " |bc )
s11=$( echo "scale=30; ( $c11 + $c12 )  /   ( $c11 * $c11 + $c11 * $c12 - 2 * $c12 * $c12 )" |bc)
s12=$( echo "scale=30; ( - $c12 )       /   ( $c11 * $c11 + $c11 * $c12 - 2 * $c12 * $c12 ) " |bc )
s44=$( echo "scale=30; ( 1 / $c44 )" | bc)
v100=$( echo "scale=30; - ( $s12/$s11 ) " |bc)
v110=$( echo "scale=30; - ( $s11 + 3 * $s12 - $s44 / 2 ) / ( 2 * $s11 + 2 * $s12 + $s44 ) " | bc)
v111=$( echo "scale=30; - ( $s11 + 2 * $s12 - $s44 / 2 ) / ( $s11 + 2 * $s12 + $s44 ) " | bc)
 
E100=$( echo "scale=30; 1/( $s11- 2 * ($s11 - $s12 - $s44 * 0.5 ) * (0) )  " |bc)
E110=$( echo "scale=30; 1/( $s11-2 * ( $s11 - $s12 - $s44 * 0.5 ) * (0.707106781^2 * 0.707106781^2 + 0.707106781^2 * 0^2 + 0^2 * 0.707106781^2 ) )   " |bc)
E111=$( echo "scale=30; 1/( $s11-2 * ( $s11 - $s12 - $s44 * 0.5 ) * (0.577350269^2 * 0.577350269^2 + 0.577350269^2 * 0.577350269^2 + 0.577350269^2 * 0.577350269^2) )   " |bc)
 
A_z=$( echo "scale=30; 2* $c44 / ($c11 - $c12)  " |bc)
A_c=$( echo "scale=30; 3* ($A_z - 1)^2 / (3* ($A_z - 1)^2 +25 * $A_z) " |bc)
A_u=$( echo "scale=30; 10* ($A_c / (1- $A_c) )  " |bc)
E_E=$( echo "scale=30; $E111/$E100  " |bc)
 
s_mo_v=$( echo "scale=30; 1 / 5 * ( $c11 - $c12 + 3 * $c44 ) " | bc)
s_mo_r=$( echo "scale=30; ( 5 * $c44 * ( $c11 - $c12 ) ) / ( 4 * $c44 + 3 * ( $c11 - $c12 ) ) " | bc)
s_mo=$( echo "scale=30; ( $s_mo_v + $s_mo_r ) / 2 " | bc)
B_mo=$( echo "scale=30; 1 / 3 * ( $c11 + 2 * $c12 ) " | bc)
y_mo=$( echo "scale=30; 9 * $B_mo * $s_mo / ( 3 * $B_mo + $s_mo ) " | bc)
iso_v=$( echo "scale=30; ( 3 * $B_mo - 2 * $s_mo )/( 2 * ( 3 * $B_mo + $s_mo ) ) " | bc)
#lat_con=$(cat ./POSCAR |head -`expr 3`| tail -1 | awk '{printf "%20.9f" , $1}')
pri_lat_con=$(cat ./CONTCAR |head -`expr 3`| tail -1 | awk '{printf "%20.9f" , $2}');  lat_con=$(echo "scale=10; $pri_lat_con * 2" | bc)
if [ $1 ]
then
echo " "
else
 
 
echo "#########################################################" >./MODULUS
echo "##                                                     ##" >>./MODULUS
echo "##  This function is useful for FCC-PRIMITIVE unitcell ##" >>./MODULUS
echo "##  Please use this function with CONTCAR and OUTCAR    ##" >>./MODULUS
echo "##  calculated by INCAR including [IBRION=6]           ##" >>./MODULUS
echo "##                                                     ##" >>./MODULUS
echo "##  Made by Ji-Hwan, Lee                               ##" >>./MODULUS
echo "##  Updated  Date: 2014.11. 07                          ##" >>./MODULUS
echo "##  Modified Date: 2016. 2. 25                          ##" >>./MODULUS
echo "#########################################################" >>./MODULUS
echo >>./MODULUS
echo >>./MODULUS
echo "#########################################################" >>./MODULUS
echo "##                  Stiffness result                   ##" >>./MODULUS
echo "#########################################################" >>./MODULUS
cat stiffness_report >> ./MODULUS
echo >>./MODULUS
echo >>./MODULUS
echo "#########################################################" >>./MODULUS
echo "##               Mechanical properties                 ##" >>./MODULUS
echo "#########################################################" >>./MODULUS
echo >>./MODULUS
echo "a0 =     " $(echo $lat_con |cut -c 1-8) >> ./MODULUS
echo "==========Compliance-Unit-(GPa)============ " >> ./MODULUS
echo "C11 =     " $(echo $c11 |cut -c 1-8) >> ./MODULUS
echo "C12 =     " $(echo $c12 |cut -c 1-8) >> ./MODULUS
echo "C44 =     " $(echo $c44 |cut -c 1-8) >> ./MODULUS
echo "==========Stiffness-Unit-(-TPa)-============ " >> ./MODULUS
echo "S11 =     " $(echo $s11 |cut -c 1-8) >> ./MODULUS
echo "S12 =     " $(echo $s12 |cut -c 1-8) >> ./MODULUS
echo "S44 =     " $(echo $s44 |cut -c 1-8) >> ./MODULUS
echo "======RESULT=====anisotropic-properties====== " >> ./MODULUS
echo "V[100] =  " $(echo $v100|cut -c 1-8)  >> ./MODULUS
echo "V[110] =  " $(echo $v110|cut -c 1-8)  >> ./MODULUS
echo "V[111] =  " $(echo $v111|cut -c 1-8)  >> ./MODULUS
echo "A_chung = " $(echo $A_c|cut -c 1-8)  >> ./MODULUS
echo "A_zener = " $(echo $A_z|cut -c 1-8)  >> ./MODULUS
echo "A_univ = " $(echo $A_u|cut -c 1-8)  >> ./MODULUS
echo "E[100]    = " $(echo $E100 |cut -c 1-4) >> ./MODULUS
echo "E[110]    = " $(echo $E110 |cut -c 1-4) >> ./MODULUS
echo "E[111]    = " $(echo $E111 |cut -c 1-4) >> ./MODULUS
echo "E111/E100 = " $(echo $E_E|cut -c 1-8)  >> ./MODULUS
echo "======RESULT========isotropic-properties====== " >> ./MODULUS
 
echo "Modul.Bulk.v    = " $(echo $b_mo_v |cut -c 1-8) >> ./MODULUS
echo "Modul.Bulk.r    = " $(echo $b_mo_r |cut -c 1-8) >> ./MODULUS
echo "Modul.Bulk(B0)  = " $(echo $b_mo|cut -c 1-8)  >> ./MODULUS
 
echo "Modul.Young.v   = " $(echo $y_mo_v |cut -c 1-8) >> ./MODULUS
echo "Modul.Young.r   = " $(echo $y_mo_r |cut -c 1-8) >> ./MODULUS
echo "Modul.Young(E0) = " $(echo $y_mo|cut -c 1-8)  >> ./MODULUS
 
echo "Modul.Shear.v   = " $(echo $s_mo_v |cut -c 1-8) >> ./MODULUS
echo "Modul.Shear.r   = " $(echo $s_mo_r |cut -c 1-8) >> ./MODULUS
echo "Modul.Shear(G0) = " $(echo $s_mo|cut -c 1-8)  >> ./MODULUS
 
echo >> ./MODULUS
echo "isotropic v =     " $(echo $iso_v |cut -c 1-8)  >> ./MODULUS
echo >>./MODULUS
echo >>./MODULUS
echo "##########################DONE############################" >>./MODULUS
 
echo lat c11  c12  c44  E100 E110 E111 v100 v110 v111 >> ./MODULUS
echo `echo $lat_con |cut -c 1-8`  `echo $c11 |cut -c 1-8`  `echo $c12 |cut -c 1-8`  `echo $c44 |cut -c 1-8`  `echo $E100 |cut -c 1-8`  `echo $E110 |cut -c 1-8`  `echo $E111 |cut -c 1-8`  `echo $v100 |cut -c 1-8`  `echo $v110 |cut -c 1-8`  `echo $v111 |cut -c 1-8`>>MODULUS
fi
rm stiff*
}
 
 
 
function HCP_mech(){
#########################################################
##                                                     ##
##  This function is useful for HCP-PRIMITIVE unitcell ##
##  Please use this function with POSCAR and OUTCAR    ##
##  calculated by INCAR including [IBRION=6]           ##
##  Add anisotropic properties of FCC (Elastic, Factor)##
##  Made by Ji-Hwan, Lee                               ##
##  Updated Date: 2016. 2. 25                          ##
##  Ref: 
##  
#########################################################
 
c11bar=0;c12bar=0;c13bar=0;c33bar=0;c66bar=0;c44bar=0; c66=0;C1=0;s11=0;s12=0;s13=0;s33=0;s44=0;s66=0;s66ref=0;
AA=0;BB=0;GG=0;CC=0;CC1=0;CC2=0;b_mo_v=0;s_mo_v=0;y_mo_v=0; b_mo_r=0;s_mo_r=0;y_mo_r=0;
b_mo=0;s_mo=0;y_mo=0; v110=0;v010=0;v001=0; iso_v=0; a_lat_con=0;c_lat_con=0;
 
egrep -10 ELASTIC OUTCAR | tail -10 >> stiffness_report
 
c11bar=$(head -n 3 ./stiffness_report | tail -n 1 | awk '{print $2}'); c11=$(echo "scale=30; ( $c11bar / 10 ) " |bc); c22=$c11
c12bar=$(head -n 3 ./stiffness_report | tail -n 1 | awk '{print $3}'); c12=$(echo "scale=30; ( $c12bar / 10 ) " |bc)
c13bar=$(head -n 3 ./stiffness_report | tail -n 1 | awk '{print $4}'); c13=$(echo "scale=30; ( $c13bar / 10 ) " |bc); c23=$c13
c33bar=$(head -n 5 ./stiffness_report | tail -n 1 | awk '{print $4}'); c33=$(echo "scale=30; ( $c33bar / 10 ) " |bc)
c66bar=$(head -n 6 ./stiffness_report | tail -n 1 | awk '{print $5}'); c66ref=$(echo "scale=30; ( $c66bar / 10 ) " |bc)
c44bar=$(head -n 7 ./stiffness_report | tail -n 1 | awk '{print $6}'); c44=$(echo "scale=30; ( $c44bar / 10 ) " |bc); c55=$c44
 
c66=$(echo "scale=30; ( ($c11-$c12)/2) " |bc)
C1=$( echo "scale=30; ( $c33 * ( $c11 + $c12 ) - 2 * $c13 ^ 2 )" |bc)
# Ref1: J.F. Nye, Physical Properties of Crystals: Their Representation by Tensors and Matrices (Oxford University Press, New York, 1985). Ch. 8 Elasticity Fourth-rank Tensor
# Ref2: H.M. Ledbetter, Journal of Physical and Chemical Reference Data 6, 1181 (1977).
s11=$(echo "scale=30; 0.5 * ( ( $c33 / $C1 ) + 1 / ( $c11 - $c12 ) ) * 1000" |bc); s22=$s11
s12=$(echo "scale=30; 0.5 * ( ( $c33 / $C1 ) - 1 / ( $c11 - $c12 ) ) * 1000" |bc)
s13=$(echo "scale=30; ( - ( $c13 / $C1 ) * 1000  )" |bc); s23=$s13
s33=$(echo "scale=30; ( ( $c11 + $c12 ) / $C1 ) * 1000  " |bc)
s44=$(echo "scale=30; ( 1 / $c44 ) * 1000 " | bc)
s66=$(echo "scale=30; 1 / $c66 * 1000 " | bc)
s66ref=$(echo "scale=30; ( $s11 - $s12 ) * 2 " |bc)
 
 
AA=$(echo "scale=30; ( $c11 + $c22 + $c33 ) / 3" |bc )
BB=$(echo "scale=30; ( $c23 + $c13 + $c12 ) / 3" |bc )
GG=$(echo "scale=30; ( $c44 + $c55 + $c66 ) / 3" |bc )
CC=$(echo "scale=30;  ( 3 * $AA + 2 * $BB + 4 * $GG ) / 5" |bc )
CC1=$(echo "scale=30; ( $AA + 4 * $BB - 2 * $GG ) / 5" |bc )
CC2=$(echo "scale=30; ( $AA - $BB + 3 * $GG ) / 5" |bc )
# s_mo : shear modulus (G)
# b_mo : bulk modulus (K)
# y_mo : young's modulus (E)
b_mo_v=$( echo "scale=30; ( 2 * $c11 + $c33 + 4* $c13 + 2 * $c12 ) / 9 " | bc)
s_mo_v=$( echo "scale=30; ( 3.5 * $c11 + $c33 - 2 * $c13 - 2.5 * $c12 + 6 * $c44 ) / 15 " | bc)
y_mo_v=$( echo "scale=30; ( $CC2 ) * ( $CC + 2 * $CC2 ) / ( $CC1 + $CC2 )" | bc)
 
b_mo_r=$( echo "scale=30;   1 / (  2 * $s11 + 1 * $s33 + 4 * $s13 +  2 *$s12 ) * 1000 " | bc) 
s_mo_r=$( echo "scale=30;  15 / ( 14 * $s11 + 4 * $s33 - 8 * $s13 - 10 *$s12 + 6 * $s44 ) * 1000 " | bc) 
# b_mo_r=$( echo "scale=30;  1 / ( $s11+$s22+$s33 + 2* ( $s23 + $s13 + $s12) )" | bc)
# s_mo_r=$( echo "scale=30; 15 / ( 4 * ( $s11+$s22+$s33 ) - 4 * ( $s23 + $s13 + $s12) + 3 * (2* $s44+$s66) )" | bc)
y_mo_r=$( echo "scale=30; 9 * $b_mo_r * $s_mo_r / (3* $b_mo_r+$s_mo_r)" | bc)
 
b_mo=$( echo "scale=30; ( $b_mo_v + $b_mo_r ) / 2" | bc)
s_mo=$( echo "scale=30; ( $s_mo_v + $s_mo_r ) / 2" | bc)
y_mo=$( echo "scale=30; ( $y_mo_v + $y_mo_r ) / 2" | bc)
 
v110=$(echo "scale=30; - ( $s12 + $s13 ) / ( 2 * $s11 )" |bc)
v010=$(echo "scale=30; - ( $s12 + $s13 ) / ( 2 * $s11 )" | bc)
v001=$(echo "scale=30; - ( $s13 ) / ( $s33 )" | bc)
 
iso_v=$( echo "scale=30; $y_mo / ( 2 * $s_mo ) - 1 " | bc)
 
#lat_con=$(cat ./POSCAR |head -`expr 3`| tail -1 | awk '{printf "%20.9f" , $1}')
a_lat_con=$(head -3 ./CONTCAR | tail -1 | awk '{printf "%20.9f", $1}'); 
c_lat_con=$(head -5 ./CONTCAR | tail -1 | awk '{printf "%20.9f", $3}' )
 
if [ $1 ]; then
return 0
else
echo "#########################################################" >./MODULUS
echo "##                                                     ##" >>./MODULUS
echo "##  This function is useful for HCP unitcell           ##" >>./MODULUS
echo "##  Please use this function with CONTCAR and OUTCAR   ##" >>./MODULUS
echo "##  calculated by INCAR including [IBRION=6]           ##" >>./MODULUS
echo "##                                                     ##" >>./MODULUS
echo "##  Made by Ji-Hwan, Lee                               ##" >>./MODULUS
echo "##  Updated  Date: 2014.11. 07                          ##" >>./MODULUS
echo "##  Modified Date: 2016. 2. 25                          ##" >>./MODULUS
echo "#########################################################" >>./MODULUS
echo >>./MODULUS
echo >>./MODULUS
echo "#########################################################" >>./MODULUS
echo "##                  Stiffness result                   ##" >>./MODULUS
echo "#########################################################" >>./MODULUS
cat stiffness_report >> ./MODULUS
echo >>./MODULUS
echo >>./MODULUS
echo "#########################################################" >>./MODULUS
echo "##               Mechanical properties                 ##" >>./MODULUS
echo "#########################################################" >>./MODULUS
echo >>./MODULUS
echo "a0 =     " $(echo $a_lat_con |cut -c 1-8) >> ./MODULUS
echo "c0 =     " $(echo $c_lat_con |cut -c 1-8) >> ./MODULUS
echo "==========Compliance-Unit-(GPa)============ " >> ./MODULUS
echo "C11 =     " $(echo $c11 |cut -c 1-8) >> ./MODULUS
echo "C12 =     " $(echo $c12 |cut -c 1-8) >> ./MODULUS
echo "C13 =     " $(echo $c13 |cut -c 1-8) >> ./MODULUS
echo "C33 =     " $(echo $c33 |cut -c 1-8) >> ./MODULUS
echo "C44 =     " $(echo $c44 |cut -c 1-8) >> ./MODULUS
echo "C66 =     " $(echo $c66 |cut -c 1-8) "C66ref =  " $(echo $c66ref|cut -c 1-8)  >> ./MODULUS
echo "==========Stiffness-Unit-(-TPa)-============ " >> ./MODULUS
echo "S11 =     " $(echo $s11 |cut -c 1-8) >> ./MODULUS
echo "S12 =     " $(echo $s12 |cut -c 1-8) >> ./MODULUS
echo "S13 =     " $(echo $s13 |cut -c 1-8) >> ./MODULUS
echo "S33 =     " $(echo $s33 |cut -c 1-8) >> ./MODULUS
echo "S44 =     " $(echo $s44 |cut -c 1-8) >> ./MODULUS
echo "S66 =     " $(echo $s66 |cut -c 1-8) "S66ref =  " $(echo $s66ref|cut -c 1-8)  >> ./MODULUS
echo "======RESULT=====anisotropic-properties====== " >> ./MODULUS
echo "V[001] =  " $(echo $v001|cut -c 1-8)  >> ./MODULUS
echo "V[010] =  " $(echo $v010|cut -c 1-8)  >> ./MODULUS
echo "V[110] =  " $(echo $v110|cut -c 1-8)  >> ./MODULUS
echo "======RESULT========isotropic-properties====== " >> ./MODULUS
 
echo "Modul.Bulk.v    = " $(echo $b_mo_v |cut -c 1-8) >> ./MODULUS
echo "Modul.Bulk.r    = " $(echo $b_mo_r |cut -c 1-8) >> ./MODULUS
echo "Modul.Bulk(B0)  = " $(echo $b_mo|cut -c 1-8)  >> ./MODULUS
 
echo "Modul.Young.v   = " $(echo $y_mo_v |cut -c 1-8) >> ./MODULUS
echo "Modul.Young.r   = " $(echo $y_mo_r |cut -c 1-8) >> ./MODULUS
echo "Modul.Young(E0) = " $(echo $y_mo|cut -c 1-8)  >> ./MODULUS
 
echo "Modul.Shear.v   = " $(echo $s_mo_v |cut -c 1-8) >> ./MODULUS
echo "Modul.Shear.r   = " $(echo $s_mo_r |cut -c 1-8) >> ./MODULUS
echo "Modul.Shear(G0) = " $(echo $s_mo|cut -c 1-8)  >> ./MODULUS
 
echo >> ./MODULUS
echo "isotropic v =     " $(echo $iso_v |cut -c 1-8)  >> ./MODULUS
echo >>./MODULUS
echo >>./MODULUS
echo "##########################DONE############################" >>./MODULUS
echo "a_Lat.Const.    ""c_Lat.Const.     ""c11           ""c12           ""c13           ""c33           ""c44           ""c66           ""iso.Poisson.  ""Modul.Bulk    ""Modul.Shear   ""Modul.Young   "  >>./MODULUS
echo $(echo $a_lat_con |cut -c 1-8)"      "$(echo $c_lat_con |cut -c 1-8)"      "$(echo $c11 |cut -c 1-8)"      "$(echo $c12 |cut -c 1-8)"      "$(echo $c13 |cut -c 1-8)"      "$(echo $c33 |cut -c 1-8)"      "$(echo $c44 |cut -c 1-8)"      "$(echo $c66 |cut -c 1-8)"      "$(echo $iso_v |cut -c 1-4)"          "$(echo $B_modu |cut -c 1-8)"      "$(echo $s_mo |cut -c 1-8)"      "$(echo $y_mo |cut -c 1-8)"      ">>./MODULUS
fi
rm stiff*
}
