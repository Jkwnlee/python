#!/bin/sh
function poscar_head(){
    num_atom=$1;    num_rlx=$2
    scale=$( echo "scale=30 ; ( $num_atom - 1 ) / 2  " | bc )
    num_fix_lay=$(echo "scale=30 ; ( $num_atom - $num_rlx ) / 2 " | bc |cut -c 1-1 )
}
 
function poscar_lattic_from_contcar_prim(){
    material_name=$(cat ./CONTCAR |head -6 |tail -1| awk '{print $1}')
    scale_num=$(cat ./CONTCAR |head -2 |tail -1)
    lat_con_a=$(cat ./CONTCAR |head -`expr 3`| tail -1 | awk '{printf "%20.9f" , $1}')
    lat_con=$(echo "scale=30 ; $lat_con_a * $scale_num " | bc )
    if [[ $lat_con ]]; then
        lat_con_a=$(cat ./CONTCAR |head -`expr 3`| tail -1 | awk '{printf "%20.9f" , $2}')
        lat_con=$(echo "scale=30 ; 2 * $lat_con_a * $scale_num " | bc )
    fi
}
 
function poscar_slab_vector_selec(){
    ax=$1;  bx=$2;  by=$3;  cz=$4
    echo $ax"         0.0000000000         0.0000000000" >> ./POSCAR
    echo $bx"         "$by"         0.0000000000" >> ./POSCAR
    echo "0.0000000000         0.0000000000         "$cz >> ./POSCAR
    echo $material_name >> ./POSCAR;      echo $num_atom >> ./POSCAR
    echo Selective Dynamics >> ./POSCAR;  echo Cartesian >> ./POSCAR
    echo "0.0000000000         0.0000000000         0.0000000000  F F F" >> ./POSCAR  ## A layer
}
 
function scaling_making_slab(){
##  scaling_making_slab $test_a $test_b $len_x $len_y $len_z $num_fix_lay $seq_num
    test_a=$1;  test_b=$2;  len_x=$3;   len_y=$4;   len_z=$5; num_fix_lay=$6 ; seq_num=$7
    cb=0;   ca=0;   cz=0
    for a in `seq 1 1 $scale`
        do
        left=$( echo  $a%$seq_num |bc )
##Define x (a) coefficient!
        ca1=$( echo "scale=30 ; $ca + $len_x " |bc )
        if [ $(bc <<< "$ca1 < $test_a") -eq 1  ]; then
            ca2=$( echo "$ca1" |cut -c 1-12)             ## WHEN ca1  < $test
            else
            ca2=$( echo "scale=30 ; $ca1 - $test_a " |bc ) ## WHEN ca1  < $test
        fi
        ca=$( echo "$ca2" |cut -c 1-12)             ## WHEN ca1  < $test
        caa=$( echo "scale=30 ; $test_a - $ca  " |bc |cut -c 1-12)
##Define y (b) coefficient!
        cb1=$( echo "scale=30 ; $cb + $len_y" |bc )
        if [ $(bc <<< "$cb1 < $test_b") -eq 1 ]; then
            cb2=$( echo "$cb1" )             ## WHEN cb1  < $test
            else
            cb2=$( echo "scale=30 ; $cb1 - $test_b " |bc )
        fi
        cb=$( echo "$cb2" |cut -c 1-12)
        cbb=$( echo "scale=30 ; $test_b - $cb  " |bc |cut -c 1-12)
        if [ $left -eq 0 ]; then
            # echo $left
            ca=0.0000000000
            cb=0.0000000000
            cbb=$( echo "scale=30 ; $test_b - $cb  " |bc |cut -c 1-12)
            caa=$( echo "scale=30 ; $test_a - $ca  " |bc |cut -c 1-12)
        fi
##Define y (c) coefficient!
        cz1=$( echo "scale=30 ; $cz + $len_z " |bc |cut -c 1-12)
        cz=$( echo "scale=30 ; $cz1"  | bc | cut -c 1-12 )
        czz=$( echo "scale=30 ; $vec_cz - $cz" |bc |cut -c 1-12)
##Print the result!
        if [ $a -le  $num_fix_lay ] ;then
            echo $ca"         "$cb"         "$cz" F F F " >> ./POSCAR
            echo $caa"         "$cbb"         "$czz" F F F" >> ./POSCAR
            else
            echo $ca"         "$cb"         "$cz" T T T" >> ./POSCAR
            echo $caa"         "$cbb"         "$czz" T T T" >> ./POSCAR
        fi
    done
}
 
function poscar_fcc_100_slab_primi() {
    ##usage: poscar_fcc_100_slab_primi 'vaccume distance' 'number of atom'
    ##please use this function with primitive FCC cell caltulation result 'CONTCAR'
 
    vaccume=$1;     num_atom=$2;    num_rlx=$3 ; seq_num=2
    poscar_head $num_atom $num_rlx; poscar_lattic_from_contcar_prim
    echo $material_name 100 $num_atom AL > ./POSCAR; echo 1.0 >> ./POSCAR
 
    vec_ax=$( echo "scale=30 ; $lat_con / sqrt(2) " | bc |cut -c 1-12)
    vec_by=$( echo "$vec_ax")
    vec_cz=$( echo "scale=30 ; $lat_con / 2 * ( $num_atom - 1 ) + $vaccume " |bc |cut -c 1-12)
    poscar_slab_vector_selec $vec_ax 0.0000000000 $vec_by $vec_cz
 
    test_a=$( echo "scale=30 ; $vec_ax " |bc )
    test_b=$( echo "scale=30 ; $vec_by " |bc )
    len_x=$( echo "scale=30 ; ($vec_ax) / 2" |bc )
    len_y=$( echo "scale=30 ; $vec_by / 2" |bc )
    len_z=$( echo "scale=30 ; $lat_con / 2 " |bc )
    cb=0;   ca=0;   cz=0
 
    scaling_making_slab $test_a $test_b $len_x $len_y $len_z $num_fix_lay $seq_num
    echo >> ./POSCAR
}
 
function poscar_fcc_110_slab_primi() {
    ##usage: poscar_fcc_110_slab_primi 'vaccume distance'  'number of atom'
        ##please use this function with primitive FCC cell caltulation result 'CONTCAR'
 
    vaccume=$1;     num_atom=$2;    num_rlx=$3;  seq_num=2
    poscar_head $num_atom $num_rlx; poscar_lattic_from_contcar_prim
 
    echo $material_name 110 $num_atom AL > ./POSCAR; echo 1.0 >> ./POSCAR
    vec_ax=$( echo "$lat_con" |cut -c 1-12)
    vec_by=$( echo "scale=30 ; $lat_con / sqrt(2) " | bc |cut -c 1-12)
    vec_cz=$( echo "scale=30 ; $lat_con / sqrt( 2 ) / 2 * ($num_atom - 1) + $vaccume " |bc |cut -c 1-12)
    poscar_slab_vector_selec $vec_ax 0.0000000000 $vec_by $vec_cz
 
    test_a=$( echo "scale=30 ; $vec_ax " |bc )
    test_b=$( echo "scale=30 ; $vec_by " |bc )
    len_x=$( echo "scale=30 ; ($vec_ax) / 2" |bc )
    len_y=$( echo "scale=30 ; $vec_by / 2" |bc )
    len_z=$( echo "scale=30 ; $lat_con  / sqrt( 2 ) / 2  " |bc )
    cb=0;   ca=0;   cz=0
 
    scaling_making_slab $test_a $test_b $len_x $len_y $len_z $num_fix_lay $seq_num
    echo >> ./POSCAR
 
}
 
function poscar_fcc_111_slab_primi(){
        ##usage: poscar_fcc_111_slab_primi 'vaccume distance' 'number of atom'
        ##please use this function with primitive FCC cell caltulation result 'CONTCAR'
 
    vaccume=$1;     num_atom=$2;    num_rlx=$3;  seq_num=3
    poscar_head $num_atom $num_rlx; poscar_lattic_from_contcar_prim
    echo $material_name 111 $num_atom AL > ./POSCAR; echo 1.0 >> ./POSCAR
 
    vec_ax=$( echo "scale=30 ; $lat_con / sqrt(2) " | bc |cut -c 1-12)
    vec_bx=$( echo "scale=30 ; $vec_ax * 0.5 " | bc |cut -c 1-12)
    vec_by=$( echo "scale=30 ; $vec_ax * 0.5 * sqrt(3) " | bc | cut -c 1-12)
    vec_cz=$( echo "scale=30 ; $lat_con * sqrt(3) / 3 * ($num_atom - 1) + $1 " |bc |cut -c 1-12)
    poscar_slab_vector_selec $vec_ax $vec_bx $vec_by $vec_cz
 
    test_a=$( echo "scale=30 ; $vec_bx + $vec_ax " |bc )
    test_b=$( echo "scale=30 ; $vec_by " |bc )
    len_x=$( echo "scale=30 ; ($vec_bx + $vec_ax) / 3" |bc )
    len_y=$( echo "scale=30 ; $vec_by / 3" |bc )
    len_z=$( echo "scale=30 ; $lat_con * sqrt(3) / 3 " |bc )
    cb=0;   ca=0;   cz=0
 
    scaling_making_slab $test_a $test_b $len_x $len_y $len_z $num_fix_lay $seq_num
    echo >> ./POSCAR
}
 
function poscar_fcc_210_slab_primi(){
        ##usage: poscar_fcc_111_slab_primi 'vaccume distance' 'number of atom'
        ##please use this function with primitive FCC cell caltulation result 'CONTCAR'
    vaccume=$1;     num_atom=$2;    num_rlx=$3 ; seq_num=10
    poscar_head $num_atom $num_rlx; poscar_lattic_from_contcar_prim
 
    echo $material_name 210 $num_atom AL > ./POSCAR; echo 1.0 >> ./POSCAR
 
    vec_ax=$( echo "scale=30 ; $lat_con * sqrt(6) * 0.5 " | bc |cut -c 1-12) #done
    vec_bx=$( echo "scale=30 ; $vec_ax * 2 / 3 " | bc |cut -c 1-12)
    vec_by=$( echo "scale=30 ; $vec_ax * sqrt(5) / 3 " | bc | cut -c 1-12)
    vec_cz=$( echo "scale=30 ; $lat_con * sqrt(5) / 10 * ($num_atom - 1) + $vaccume " |bc |cut -c 1-12)
    poscar_slab_vector_selec $vec_ax $vec_bx $vec_by $vec_cz
 
    test_a=$( echo "scale=30 ; $vec_bx + $vec_ax " |bc )
    test_b=$( echo "scale=30 ; $vec_by " |bc )
    len_x=$( echo "scale=30 ; ($vec_bx + $vec_ax ) * 0.7" |bc )
    len_y=$( echo "scale=30 ; ($vec_by ) * 0.7" |bc )
    len_z=$( echo "scale=30 ; $lat_con * sqrt(5) / 10 " |bc )
    cb=0;   ca=0;   cz=0
 
    scaling_making_slab $test_a $test_b $len_x $len_y $len_z $num_fix_lay $seq_num
 
    echo >> ./POSCAR
}
 
 
 
function poscar_fcc_211_slab_primi(){
        ##usage: poscar_fcc_111_slab_primi 'vaccume distance' 'number of atom'
        ##please use this function with primitive FCC cell caltulation result 'CONTCAR'
 
    vaccume=$1;     num_atom=$2;    num_rlx=$3;  seq_num=6
    poscar_head $num_atom $num_rlx; poscar_lattic_from_contcar_prim
    echo $material_name 211 $num_atom AL > ./POSCAR; echo 1.0 >> ./POSCAR
 
    vec_ax=$( echo "scale=30 ; $lat_con * sqrt(2) * 0.5 " | bc |cut -c 1-12) #done
    vec_by=$( echo "scale=30 ; $lat_con * sqrt(3) " | bc | cut -c 1-12)
    vec_cz=$( echo "scale=30 ; $lat_con * sqrt(6) * 0.5 / 6 * ($num_atom - 1) + $vaccume " |bc |cut -c 1-12)
    poscar_slab_vector_selec $vec_ax 0.0000000000 $vec_by $vec_cz
 
    test_a=$( echo "scale=30 ; $vec_ax " |bc )
    test_b=$( echo "scale=30 ; $vec_by " |bc )
    len_x=$( echo "scale=30 ; ($vec_ax ) * 0.5 " |bc )
    len_y=$( echo "scale=30 ; ($vec_by ) * 2 / 3" |bc )
    len_z=$( echo "scale=30 ; $lat_con * sqrt(6) / 2 / 6  " |bc )
    cb=0;   ca=0;   cz=0
 
    scaling_making_slab $test_a $test_b $len_x $len_y $len_z $num_fix_lay $seq_num
    echo >> ./POSCAR
 
}
 
 
function poscar_fcc_311_slab_primi(){
        ##usage: poscar_fcc_111_slab_primi 'vaccume distance' 'number of atom'
        ##please use this function with primitive FCC cell caltulation result 'CONTCAR'
 
    vaccume=$1;     num_atom=$2;    num_rlx=$3;  seq_num=11
    poscar_head $num_atom $num_rlx; poscar_lattic_from_contcar_prim
    echo $material_name 311 $num_atom AL > ./POSCAR; echo 1.0 >> ./POSCAR
 
    vec_ax=$( echo "scale=30 ; $lat_con * sqrt(6) / 2 " | bc |cut -c 1-12) #done
    vec_bx=$( echo "scale=30 ; $lat_con * sqrt(6) / 2 * 5 / 6 " | bc |cut -c 1-12)  ## Cos (theta) =5/6
    vec_by=$( echo "scale=30 ; $lat_con * sqrt(6) / 2 * sqrt(11) / 6  " | bc | cut -c 1-12)
    vec_cz=$( echo "scale=30 ; $lat_con / sqrt(11)  * ($num_atom - 1) + $vaccume " |bc |cut -c 1-12)
    poscar_slab_vector_selec $vec_ax $vec_bx $vec_by $vec_cz
 
    test_a=$( echo "scale=30 ; $vec_bx + $vec_ax " |bc )
    test_b=$( echo "scale=30 ; $vec_by " |bc )
    len_x=$( echo "scale=30 ; ($vec_ax + $vec_bx) *  0.2727 " |bc )
    len_y=$( echo "scale=30 ; $vec_by * 0.2727 " |bc )
    len_z=$( echo "scale=30 ; $lat_con / sqrt(11) " |bc )
    cb=0;   ca=0;   cz=0
 
    scaling_making_slab $test_a $test_b $len_x $len_y $len_z $num_fix_lay $seq_num
 
    echo >> ./POSCAR
}
 
 
function poscar_fcc_331_slab_primi(){
        ##usage: poscar_fcc_111_slab_primi 'vaccume distance' 'number of atom'
        ##please use this function with primitive FCC cell caltulation result 'CONTCAR'
 
    vaccume=$1;     num_atom=$2;    num_rlx=$3;  seq_num=19 #seq_num=11
    poscar_head $num_atom $num_rlx; poscar_lattic_from_contcar_prim
    echo $material_name 311 $num_atom AL > ./POSCAR; echo 1.0 >> ./POSCAR
 
    vec_ax=$( echo "scale=30 ; $lat_con * sqrt(10) / 2 " | bc |cut -c 1-12) #done
    vec_bx=$( echo "scale=30 ; $lat_con * sqrt(10) / 2 * 0.9 " | bc |cut -c 1-12)
    vec_by=$( echo "scale=30 ; $lat_con * sqrt(10) / 2 * sqrt(0.19) " | bc | cut -c 1-12)
    vec_cz=$( echo "scale=30 ; $lat_con / sqrt(19)  * ($num_atom - 1) + $vaccume " |bc |cut -c 1-12)
    poscar_slab_vector_selec $vec_ax $vec_bx $vec_by $vec_cz
 
    test_a=$( echo "scale=30 ; $vec_bx + $vec_ax " |bc )
    test_b=$( echo "scale=30 ; $vec_by " |bc )
    len_x=$( echo "scale=30 ; ($vec_ax + $vec_bx) * ( 1 - sqrt(10) / 10 ) " |bc )
    len_y=$( echo "scale=30 ; $vec_by * ( 1 - sqrt(10) / 10 )  " |bc )
    len_z=$( echo "scale=30 ; $lat_con / sqrt(19) " |bc )
 
    scaling_making_slab $test_a $test_b $len_x $len_y $len_z $num_fix_lay $seq_num
    echo >> ./POSCAR
}
 
 
 
function atom_poscar(){
cat >POSCAR<<qwer
$1 - atom
1
15 0 0
0 15.5 0
0 0 16
$1
1
Direct
0 0 0
qwer
 
}
 
function atom_kpoint(){
cat >KPOINTS<<qwer
Automatic mesh
0
Gamma
  1.0     1.0     1.0
  0.25    0.25    0.25
 
qwer
}
 
function fcc_poscar_primitive_ev(){
cat > ./POSCAR<<qwer
$1 - SC
$3
0 $2 $2
$2 0 $2
$2 $2 0
$1
1
Direct
0.0000000000000000  0.0000000000000000  0.0000000000000000
qwer
 
}
 
 
function fcc_poscar_primitive() {
cat > ./POSCAR <<qwer
FCC Primitive cell $name 
1.00000000000000     
0 $2 $2
$2 0 $2
$2 $2 0
$1
1
Direct
0.0000000000000000  0.0000000000000000  0.0000000000000000
 
qwer
}
 
 
function fcc_poscar_100_base() {
cat > 100POSCAR <<qwer
$1 FCC 100-oriented
1.0
$2 0.0000000000 0.0000000000
0.0000000000 $3 0.0000000000
0.0000000000 0.0000000000 $4
$1
4
Direct
0 0.5 0.5
0.5 0 0.5
0.5 0.5 0
0 0 0
 
qwer
}
 
function fcc_poscar_110_base() {
cat > 110POSCAR <<qwer
$1 FCC 110-oriented
1.0
$2 0.0000000000 0.0000000000
0.0000000000 $3 0.0000000000
0.0000000000 0.0000000000 $4
$1
2
Direct
0.5 0.5 0.5
0 0 0
qwer
 
}
 
function fcc_poscar_111_base() {
cat > 111POSCAR <<qwer
$1 FCC 111-oriented
1.0
$2         0.0000000000         0.0000000000
$3         $4         0.0000000000
0.0000000000         0.0000000000        $5
$1
3
Direct
0.666700006         0.666700006         0.333299994
0.333299994         0.333299994         0.666700006
0.000000000         0.000000000         0.000000000 
 
qwer
}
 
 
function fcc_poscar_primitive_ev(){
cat > ./POSCAR <<qwer
$1 energy- volume graph
$3
0 $2 $2
$2 0 $2
$2 $2 0
$1
1
Direct
0.0000000000000000  0.0000000000000000  0.0000000000000000
 
qwer
}
 
 
function dia_poscar_primitive() { 
cat > ./POSCAR <<qwer
$1
1.00000000000000
$2    $2    0.0000000000000000
0.0000000000000000    $2    $2
$2    0.0000000000000000    $2
$1
2
Direct
0.0000000000000000  0.0000000000000000  0.0000000000000000
0.2500000000000000  0.2500000000000000  0.2500000000000000
 
qwer
}
 
 
function dia_poscar_111_base() {
cat > 111POSCAR <<qwer
$1 111-oriented Diamond structure
1.0
$2         0.0000000000         0.0000000000
$3         $4         0.0000000000
0.0000000000         0.0000000000         $5
$1
6
Direct
0.000000000         0.000000000         0.250000000
0.666700006         0.666700006         0.333299994
0.666700006         0.666700006         0.583299994
0.333299994         0.333299994         0.666700006
0.333299994         0.333299994         0.916700006
0.000000000         0.000000000         0.000000000
 
qwer
}
 
function dia_poscar_100_base() {
cat > 100POSCAR <<qwer
$1 100-oriented Diamond structure
1.0
$2        0.0000000000         0.0000000000
0.0000000000        $3         0.0000000000
0.0000000000         0.0000000000        $4
$1
8
Direct
0.250000000         0.250000000         0.250000000
0.000000000         0.500000000         0.500000000
0.500000000         0.000000000         0.500000000
0.750000000         0.250000000         0.750000000
0.250000000         0.750000000         0.750000000
0.000000000         0.000000000         0.000000000
0.500000000         0.500000000         0.000000000
0.750000000         0.750000000         0.250000000
 
qwer
}
 
function dia_poscar_110_base() {
cat > 110POSCAR <<qwer
$1 110-oriented Diamond structure
1.0
$2 0.0000000000 0.0000000000
0.0000000000 $3 0.0000000000
0.0000000000 0.0000000000 $4
$1
4
Direct
0.000000000         0.000000000         0.000000000
0.500000000         0.500000000         0.500000000
0.750000000         0.000000000         0.500000000
0.250000000         0.500000000         0.000000000
 
qwer
}
 
 
function hcp_poscar(){
#if [[ "$1" == `  ` ]]
#then
#echo Warning Follow the following comment,  >> POSCAR
#echo "hcp_poscar \"name_of_material\" \"a_axis_length\" \"c_axis_length\" " >> POSCAR
#fi
 
a_axis_leng=$( echo " $2 " | bc -l |awk '{printf "%20.15f", $1}' )
b_x_axis_leng=$( echo " - $a_axis_leng / 2 " | bc -l | awk '{printf "%20.15f", $1}' )
b_y_axis_leng=$( echo " $a_axis_leng * sqrt( 3 ) / 2 " | bc -l|awk '{printf "%20.15f", $1}' )
cat > ./POSCAR <<qwer
HCP Unit cell of $1
1.000000000000000
$a_axis_leng   0.000000000000000   0.000000000000000
$b_x_axis_leng   $b_y_axis_leng   0.000000000000000
0.000000000000000   0.000000000000000   $3
$1
2
Direct
0.333333333333333  0.666666666666667  0.250000000000000
0.666666666666667  0.333333333333333  0.750000000000000
 
qwer
 
}
function bcc_poscar_primitive() {
cat > ./POSCAR <<qwer
BCC Primitive cell of $1
1.00000000000000
-$2  $2  $2
 $2 -$2  $2
 $2  $2 -$2
$1
1
Direct
0.0000000000000000  0.0000000000000000  0.0000000000000000
 
qwer
}
