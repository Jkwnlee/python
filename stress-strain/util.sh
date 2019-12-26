function Vol_Cal(){
if [[ $1 == 1 ]]
  then
  sys="CONTCAR"
  x=3
  unit="Aungstrom^3"
  scale=$(head -`expr $x - 1` $sys |tail -1 | awk '{printf $1}')
else
  echo "If the basis-set of calculation is VASP, type "v", or if one of calculation is Dmol3, type "d" "\c
  read name1
  if [ "$name1" = "v" ]
     then
     sys="POSCAR"
     x=3
     unit="Aungstrom^3"
     scale=$(head -`expr $x - 1` $sys |tail -1 | awk '{printf $1}')
     echo "start to calculate Volume of your "$sys
  elif [ "$name1" = "d" ]
    then
    sys="INCOOR"
    x=2
    unit="Bohr^3"
    scale=1
  fi
fi
 
ax=$(head -`expr $x` $sys |tail -1 | awk '{printf $1}')
ay=$(head -`expr $x` $sys |tail -1 | awk '{printf $2}')
az=$(head -`expr $x` $sys |tail -1 | awk '{printf $3}')
bx=$(head -`expr $x + 1` $sys |tail -1 | awk '{printf $1}')
by=$(head -`expr $x + 1` $sys |tail -1 | awk '{printf $2}')
bz=$(head -`expr $x + 1` $sys |tail -1 | awk '{printf $3}')
cx=$(head -`expr $x + 2` $sys |tail -1 | awk '{printf $1}')
cy=$(head -`expr $x + 2` $sys |tail -1 | awk '{printf $2}')
cz=$(head -`expr $x + 2` $sys |tail -1 | awk '{printf $3}')
 
Volume=$(echo "scale=30; (($ay * $bz - $az * $by ) * $cx + ($az * $bx - $ax * $bz ) * $cy + ($ax * $by - $ay * $bx ) * $cz ) * $scale * $scale * $scale " | bc | awk '{printf "%22.19f",$1}')
 
}
 
 
 
 
 
 
 
#!/bin/sh
# Written date: 2015.04.05 by Ji-Hwan Lee (jkwnlee88@gmail.com)
# Edited date : 2015.08.04 by Ji-Hwan Lee (jkwnlee88@gmail.com)
# Purpose :   Rearrange nodes in KOHN for using
# function  : Manually change the information of nodes in script
#                       Time following the introduction
# Usage   : type "rewrite_job"
#example: withing 
function rewrite_job(){
if [[ `head -1 job*sh` == 0 ]]; then
echo "there is no jobscript to edit"
return 0
fi
if [[ $1 -eq ' ' ]]; then
        num_lin=`wc -l job*.sh |awk '{print $1-2}'`
        tail -$num_lin job*.sh > re_job_proc
        head -1 job*.sh > new_job.sh
        name1=99
        echo "Choose number of node and core, ex: for x004:ppn12, type 4 12 "
        read name1 name2
        while [[ $name1 ]]; do
           selec_node $name1 $name2 ; echo $a1 >> re_job_proc1
           echo "Do u need more node and core? if yes please type node & core, if not please type ENTER"
           read name1 name2
        done
        nodes=`transpose_job re_job_proc1`
        echo $head"="$nodes >> new_job.sh
        walltm=`grep "walltime="  re_job_proc`
        test_par=`echo $walltm |cut -c 25-26`
        if [[ $test_par -eq '0' ]]; then
           sed -i "s/$walltm/#PBS -l walltime=$time:00:00/g" re_job_proc
        else
           echo "#PBS -l walltime=$time:00:00" >> new_job.sh
        #   echo $a
        fi
        cat re_job_proc >> new_job.sh
        rm re_job_proc re_job_proc1 re_job_proc2
else
        num_lin=`wc -l job*.sh |awk '{print $1-2}'`
        head -1 job*.sh > new_job.sh
        tail -$num_lin job*.sh > re_job_proc
        selec_node $1
        echo $head"="$a1 >> new_job.sh
        cat re_job_proc  >> new_job.sh
        rm re_job_proc re_job_proc1
fi
}
 
function transpose_job()
{
 
declare -a array=( )                      # we build a 1-D-array
 
read -a line < "$1"                       # read the headline
 
COLS=${#line[@]}                          # save number of columns
 
index=0
while read -a line; do
    for (( COUNTER=0; COUNTER<${#line[@]}; COUNTER++ )); do
        array[$index]=${line[$COUNTER]}
        ((index++))
    done
done < "$1"
 
for (( ROW = 0; ROW < COLS; ROW++ )); do
  for (( COUNTER = ROW; COUNTER < ${#array[@]}; COUNTER += COLS )); do
    printf "%s" ${array[$COUNTER]}
    if [ $COUNTER -lt $(( ${#array[@]} - $COLS )) ]
    then
        printf "+"
    fi
  done
  printf "\n"
done
}
 
 
function selec_node(){
## Call values for different nodes, core, and walltime in "$1", "core", and "time"
  function core_cond(){
    if [[ $2 -le 0 ]]; then
      if [[ $1 -le 8 ]]; then
        core=8
      else
        core=12
      fi
    else
      core=`echo $2`
    fi
  }
 
time=999
head=`echo \#PBS -l nodes`
if [[ $1 == 0 ]]; then
  print "There is no input, please type "selec_node \$number_of_node \$number_of_cores" "
else
        if [[ $1 ==    1 ]]; then
                core_cond $1 $2
                a1=`echo x001\:ppn=$core`
                time=24
        elif [[ $1 ==  2 ]]; then
                core_cond $1 $2
                a1=`echo x002\:ppn=$core`
                time=24
        elif [[ $1 ==  3 ]]; then
                core_cond $1 $2
                a1=`echo x003\:ppn=$core`
                time=24
        elif [[ $1 ==  4 ]]; then
                core_cond $1 $2
                a1=`echo x004\:ppn=$core`
                time=24
        elif [[ $1 ==  5 ]]; then
                core_cond $1 $2
                a1=`echo x005\:ppn=$core`
                time=48
        elif [[ $1 ==  6 ]]; then
                core_cond $1 $2
                a1=`echo x006\:ppn=$core`
                time=48
        elif [[ $1 ==  7 ]]; then
                core_cond $1 $2
                a1=`echo x007\:ppn=$core`
                time=48
        elif [[ $1 ==  8 ]]; then
                core_cond $1 $2
                a1=`echo x008\:ppn=$core`
                time=48
        elif [[ $1 ==  9 ]]; then
                core_cond $1 $2
                a1=`echo x009\:ppn=$core`
                time=72
        elif [[ $1 == 10 ]]; then
                core_cond $1 $2
                a1=`echo x010\:ppn=$core`
                time=48
        elif [[ $1 == 11 ]]; then
                core_cond $1 $2
                a1=`echo x011\:ppn=$core`
                time=48
        elif [[ $1 == 12 ]]; then
                core_cond $1 $2
                a1=`echo x012\:ppn=$core`
                time=48
        elif [[ $1 == 13 ]]; then
                core_cond $1 $2
                a1=`echo x013\:ppn=$core`
                time=999
        elif [[ $1 == 14 ]]; then
                core_cond $1 $2
                a1=`echo x014\:ppn=$core`
                time=999
        elif [[ $1 == 15 ]]; then
                core_cond $1 $2
                a1=`echo x015\:ppn=$core`
                time=999
        elif [[ $1 == 16 ]]; then
                core_cond $1 $2
                a1=`echo x016\:ppn=$core`
                time=999
        elif [[ $1 == 17 ]]; then
                core_cond $1 $2
                a1=`echo x017\:ppn=$core`
                time=999
        elif [[ $1 == 18 ]]; then
                core_cond $1 $2
                a1=`echo x018\:ppn=$core`
                time=999
        elif [[ $1 == 19 ]]; then
                core_cond $1 $2
                a1=`echo x019\:ppn=$core`
                time=999
        fi
fi
}
 
 
 
 
 
 
 
 
 
function poscar_direc2cart (){
  if [[ `ls -l |grep CONTCAR | head -1 | awk '{print $9}'` == "CONTCAR" ]]; then
     filename="CONTCAR"
    echo "CONTCAR is in this folder"
    echo "Let's name cartesian type CONTCAR with name - new_CONTCAR"
  else
    echo "In this folder, there is no CONTCAR"
    echo "please type the name of POSCAR type filename"
    read filename   
    echo "Let's name cartesian type "$filename" with name - new_"$filename
  fi
 
  declare -a len
  # filename="CONTCAR"
  IFS=$'\n' xyz=($(cat "$filename" ))
  IFS=$' ' x=($(echo ${xyz[2]}))
  IFS=$' ' y=($(echo ${xyz[3]}))
  IFS=$' ' z=($(echo ${xyz[4]}))
  for (( i = 0; i < ${#x[@]} ; i++ )); do
    len["$i"]=`echo  ${x["$i"]} + ${y["$i"]} + ${z["$i"]}   |bc`
    # echo ${len["$i"]}
  done
  lin_num=$(wc -l "$filename" |awk '{print $1}') 
  num_atom=$(cat "$filename" | head -7 |tail -1|  awk '{print $1 + $2 + $3 + $4 + $5 + $6 }' )
  if [[ `echo ${xyz[7]} | cut -c 1-1` == "S" ]]; then
    # echo ${xyz[7]}
    head -8 "$filename" > new_$filename
                echo "Cartesian" >> new_$filename
    head -`expr $num_atom + 9` "$filename" | tail -`expr $num_atom` | awk '{printf "  %18.16f  ",$1*'${len[0]}'}  {printf "%18.16f  ",$2*'${len[1]}'} {printf "%18.16f   ",$3*'${len[2]}'} {print $4"  "$5"  "$6}' >> new_$filename
  else
    head -8 "$filename" > new_$filename
                echo "Cartesian" >> new_$filename
    head -`expr $num_atom + 8` "$filename" | tail -`expr $num_atom` | awk '{printf "  %18.16f  ",$1*'${len[0]}'}  {printf "%18.16f  ",$2*'${len[1]}'} {printf "%18.16f   ",$3*'${len[2]}'} {print $4"  "$5"  "$6}' >> new_$filename
    # echo it is not selective 
  fi
 
}
 
 
 
 
function POTCAR_PBE_make (){
## Making       : Jihwan Lee
## Date         : 2015. 07. 23
## Purpose      : Making 1st generated  POTCAR for a element with a condition which require the biggest number of considered electron in POTCAR_list
 
path_pot=`echo /opt/VASP/POTCAR/PAW_PBE`
path_pot_list=`echo /home/jihwan/vasp/vasp_utility/jihwan_script/POTCAR_list`
 
name4=`grep "  $4  " $path_pot_list |head -1 |cut -c 9-15 `
name3=`grep "  $3  " $path_pot_list |head -1 |cut -c 9-15`
name2=`grep "  $2  " $path_pot_list |head -1 |cut -c 9-15`
name1=`grep "  $1  " $path_pot_list |head -1 |cut -c 9-15`
if [[ $4 == "" ]]; then
  if [[ $3 == "" ]]; then
    if [[ $2 == "" ]]; then
      if [[ $1 == "" ]]; then
        echo errrrrr
      else
        cat "$path_pot"/"$name1"/POTCAR > ./POTCAR
      fi
    else
      cat "$path_pot"/"$name1"/POTCAR "$path_pot"/"$name2"/POTCAR > ./POTCAR
    fi
  else
    cat "$path_pot"/"$name1"/POTCAR "$path_pot"/"$name2"/POTCAR "$path_pot"/"$name3"/POTCAR > ./POTCAR
  fi 
else
  cat "$path_pot"/"$name1"/POTCAR "$path_pot"/"$name2"/POTCAR "$path_pot"/"$name3"/POTCAR "$path_pot"/"$name4"/POTCAR > ./POTCAR
fi
 
grep TIT POTCAR
grep ZVAL POTCAR
grep ENMAX POTCAR
}
 
 
function POTCAR_LDA_make (){
## Making       : Jihwan Lee
## Date         : 2015. 07. 23
## Purpose      : Making 1st generated  POTCAR for a element with a condition which require the biggest number of considered electron in POTCAR_list
 
path_pot=`echo /opt/VASP/POTCAR/PAW_LDA`
path_pot_list=`echo /home/jihwan/vasp/vasp_utility/jihwan_script/POTCAR_list`
 
name4=`grep "  $4  " $path_pot_list |head -1 |cut -c 9-15 `
name3=`grep "  $3  " $path_pot_list |head -1 |cut -c 9-15`
name2=`grep "  $2  " $path_pot_list |head -1 |cut -c 9-15`
name1=`grep "  $1  " $path_pot_list |head -1 |cut -c 9-15`
if [[ $4 == "" ]]; then
  if [[ $3 == "" ]]; then
    if [[ $2 == "" ]]; then
      if [[ $1 == "" ]]; then
        echo errrrrr
      else
        cat "$path_pot"/"$name1"/POTCAR > ./POTCAR
      fi
    else
      cat "$path_pot"/"$name1"/POTCAR "$path_pot"/"$name2"/POTCAR > ./POTCAR
    fi
  else
    cat "$path_pot"/"$name1"/POTCAR "$path_pot"/"$name2"/POTCAR "$path_pot"/"$name3"/POTCAR > ./POTCAR
  fi 
else
  cat "$path_pot"/"$name1"/POTCAR "$path_pot"/"$name2"/POTCAR "$path_pot"/"$name3"/POTCAR "$path_pot"/"$name4"/POTCAR > ./POTCAR
fi
 
grep TIT POTCAR
grep ZVAL POTCAR
grep ENMAX POTCAR
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
 
 
 
function gen_kpoint(){
  # gen_kpoint F 12 12 1
  if [[ $1 == T ]]; then
    Gamma="Gamma"
  else
    Gamma="Monkhorst-Pack"
  fi
  cat >KPOINTS<<qwer
Automatic mesh
0
$Gamma
$2 $3 $4
0.0 0.0 0.0
 
qwer
}
 
 
# Please start with "echo "IGOR" >> "$name".itx"
function make_itx (){
name=$1
dataset_name=$2
wave1=$3
wave2=$4
wave3=$5
wave4=$6
# echo "IGOR" >> "$name".itx
echo "WAVES/D" $wave1"    "$wave2"    "$wave3"    "$wave4 >> "$name".itx
echo "BEGIN" >> "$name".itx
cat ./"$dataset_name" >> "$name".itx
echo "END" >> "$name".itx
echo >> "$name".itx
}
    
 
 
 
 
function read_poscar(){
  ax=$(head -`expr 3` "$1"/POSCAR |tail -1 | awk '{printf $1}')
  ay=$(head -`expr 3` "$1"/POSCAR |tail -1 | awk '{printf $2}')
  az=$(head -`expr 3` "$1"/POSCAR |tail -1 | awk '{printf $3}')
  bx=$(head -`expr 4` "$1"/POSCAR |tail -1 | awk '{printf $1}')
  by=$(head -`expr 4` "$1"/POSCAR |tail -1 | awk '{printf $2}')
  bz=$(head -`expr 4` "$1"/POSCAR |tail -1 | awk '{printf $3}')
  cx=$(head -`expr 5` "$1"/POSCAR |tail -1 | awk '{printf $1}')
  cy=$(head -`expr 5` "$1"/POSCAR |tail -1 | awk '{printf $2}')
  cz=$(head -`expr 5` "$1"/POSCAR |tail -1 | awk '{printf $3}')
  len_a=$(echo "scale=30 ; sqrt( $ax ^2 + $ay ^2 + $az ^2) " | bc |cut -c 1-18) #done
  len_b=$(echo "scale=30 ; sqrt( $bx ^2 + $by ^2 + $bz ^2) " | bc |cut -c 1-18) #done
  len_c=$(echo "scale=30 ; sqrt( $cx ^2 + $cy ^2 + $cz ^2) " | bc |cut -c 1-18) #done
  echo $len_a $len_b $len_c
  if [[ `ls -al |grep OUTCAR` ]]; then
    kd_a=$(grep " generate k-points for:" $1/OUTCAR |head -1 |awk '{print $4}')
    kd_b=$(grep " generate k-points for:" $1/OUTCAR |head -1 |awk '{print $5}')
    kd_c=$(grep " generate k-points for:" $1/OUTCAR |head -1 |awk '{print $6}')
    echo $kd_a $kd_b $kd_c
    kd_aa=$(echo "scale=30 ; 1 / $len_a /$kd_a *3.14159265359*2" | bc |cut -c 1-18)
    kd_bb=$(echo "scale=30 ; 1 / $len_b /$kd_b *3.14159265359*2" | bc |cut -c 1-18)
    kd_cc=$(echo "scale=30 ; 1 / $len_b /$kd_c *3.14159265359*2" | bc |cut -c 1-18)
    echo $kd_aa $kd_bb $kd_cc
  fi
}
 
 
function kpoint_slab(){
  # $1: where the bulk poscar
  # gen_kpoint F 12 12 1
  # kpoint_slab ../1_vc_opt/ . T
 
  read_poscar $1
  bulk_kd_a=$kd_aa
  bulk_kd_b=$kd_bb
  kspacing=0.125
 
  if [[ $3 ]]; then
    read_poscar .
    kpa=$(echo "scale=30 ; 1 / $len_a / $bulk_kd_a * 3.14159265359 * 2" | bc | awk '{printf "  %18.0f  ",$1}')
    kpb=$(echo "scale=30 ; 1 / $len_b / $bulk_kd_b * 3.14159265359 * 2" | bc | awk '{printf "  %18.0f  ",$1}')
    gen_kpoint $2 $kpa $kpb 1
    sed -i "s/KSPACING/#KSPACING/g"  ./INCAR
    sed -i "s/KGAMMA/#KGAMMA/g"  ./INCAR
  else
    read_poscar $3
    kpa=$(echo "scale=30 ; 1 / $len_a / $bulk_kd_a * 3.14159265359 * 2" | bc | awk '{printf "  %18.0f  ",$1}')
    kpb=$(echo "scale=30 ; 1 / $len_b / $bulk_kd_b * 3.14159265359 * 2" | bc | awk '{printf "  %18.0f  ",$1}')
    gen_kpoint $2 $kpa $kpb 1
    sed -i "s/KSPACING/#KSPACING/g"  $2/INCAR
    sed -i "s/KGAMMA/#KGAMMA/g"  $2/INCAR
  fi
 
 
}
