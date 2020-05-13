#!bin/sh
export ISODATA=~/findsym/isobyu/
export findsym_path=~/findsym/bin

$findsym_path/findsym_cifinput $1 > $2
$findsym_path/findsym $2 > $3

