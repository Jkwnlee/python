pathutil=/home/jihwan/vasp/vasp_utility/jihwan_script/
. $pathutil/make_igor_itx.sh
function make_xy_data(){
data_name=$1
cat $data_name |awk '{print $1}' > ./Adata.txt
cat $data_name |awk '{print $2}' > ./Cdata.txt
cat $data_name |awk '{print $5}' > ./edata.txt
cat $data_name |awk '{print $1" "$2" "$5}' > ./data.txt
python $pathutil/interpolate/max_min.py > ./xy.data
 rm Adata.txt ./Cdata.txt ./edata.txt ./data.txt
}
 
function temp(){
make_xy_data 3_EV_2D/EV_DATA
 
$pathutil/interpolate/1.interpolate_program/interpolate.x
. $pathutil/interpolate/1.interpolate_program/remove_line.sh
rm interpolation-grid.data plot.data xfarbe.data xfarbe.sym.grid
# len=`wc -l wave.data | awk '{print $1}'`
 
# ii=  ##function name
potcar=$1
functional=$2
rm "$potcar"_"$functional".itx
echo "IGOR" >> "$potcar"_"$functional".itx
cat ./2_EV/c_o_a_data | awk '{print $1" \t "$2}'  > data1
make_itx "$potcar"_"$functional"  data1 EV_a_"$potcar"_"$functional"  EV_c_"$potcar"_"$functional"
cat ./plot.data2 | awk '{print $1" \t "$2" \t "$3}' > data1
make_itx "$potcar"_"$functional"  data1 EV2_a_"$potcar"_"$functional" EV2_c_"$potcar"_"$functional" EV2_e_"$potcar"_"$functional"
rm data1 plot.data2 
 rm xy.data
 
 
sed "s/potcar/$potcar/g" $pathutil/interpolate/2.contour_macro_igor/macro_igor.txt > "$potcar"_"$functional"1
sed "s/functional/$functional/g" "$potcar"_"$functional"1 > ./"$potcar"_"$functional"
for leng in `seq 1 1 $(wc -l "$potcar"_"$functional" |awk '{print $1}')` ; do
echo X `head -$leng "$potcar"_"$functional" |tail -1 ` >> "$potcar"_"$functional".itx
done
 
rm "$potcar"_"$functional"1 "$potcar"_"$functional"
 
}
