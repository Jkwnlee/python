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
