#!/bin/bash
d=`dirname $0`
#Name of the data directory
data=$1

DIR="$data/tmp" #DIR="$d/$data/tmp"
DIR_1="$data/tmp/" #DIR_1="$d/$data/tmp/"

mkdir -p $DIR
#Arguments. User can modify these arguments.

#Input data file name.
input="$data/input.txt" #input="$d/$data/input.txt"

#Input banlist file
ban="$data/banlist.txt"

#Input whitelist file
white="$data/whitelist.txt"

#Max parents. We set the number of maximum parents to 2 
maxparent=2

#k is the number of structure considered in each step of k-best structure learning algorithm. We set the value of k=100
k=100

#Model averaging threshold. We set the model averaging threshold to 0.5
THR=0.5

echo "test: \$input: $input; \$ban: $ban; \$white: $white; \$DIR: $DIR; \$DIR_1: $DIR_1"

#cd ./tools/
echo "Dataset with number of variables:"
head -1 $input|wc -w
echo "Dataset with number of samples:"
wc -l $input|cut -d' ' -f1
echo "Parameters are: max parent = $maxparent, k= $k and model averaging threshold = $THR"

start_time1=`date +%s`
#execute local score
./tools/network_score $input $ban $white $maxparent $DIR
end_time=`date +%s`
echo "Local score execution time was `expr $end_time - $start_time1` s."

start_time=`date +%s`
#execute k-best structure learning
#sh $d/k-best/src/data2netk_poster.sh $input $DIR $k $maxparent
sh ./tools/k-best/src/data2netk_poster.sh $input $DIR $k $maxparent

end_time=`date +%s`
echo "k-best parent total execution time was `expr $end_time - $start_time` s."

start_time=`date +%s`
#execute print structure
#./structure $input $DIR/postProbEachEdge.txt $d/$data/model_averaging_probabilities.txt $d/$data/model_structure.txt $THR
./tools/structure $input $DIR/postProbEachEdge.txt $data/model_averaging_probabilities.txt $data/model_structure.txt $THR

end_time1=`date +%s`
echo "Preparing structure output file from model averaging matrix. Execution time was `expr $end_time1 - $start_time` s."
echo "Total execution time was `expr $end_time1 - $start_time1` s."
rm -r $DIR_1*
rmdir $DIR
#cd-
