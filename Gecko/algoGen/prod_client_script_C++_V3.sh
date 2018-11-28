#!/bin/bash

      #1pathdata 2ppathlog 3phidneuron  4pmethod 5pgen  6pindiv     7pkmer 8shufTrainTestDataType 9noisefactor 10detailres

timestamp=$(($(date +%s%N)))
START_TIME=$SECONDS
#python clientNN.py $timestamp &
echo "nproc="
nproc
echo "OMP_NUM_THREADS = $OMP_NUM_THREADS"
if [ -z "$OMP_NUM_THREADS" ];then
    mpinpproc=$(nproc)
else
    mpinpproc=$OMP_NUM_THREADS
fi

echo "mpirun -np $mpinpproc python3 clientNNc++.py $timestamp $1"
mpirun -np $mpinpproc python3 clientNNc++.py $timestamp $1  &
#client_pid=$!
echo id = $timestamp
echo $@

timestamp1=$(awk 'BEGIN {srand(); print srand()}')
echo $timestamp $1
Producteurv2/farmer2   $timestamp $1
timestamp2=$(awk 'BEGIN {srand(); print srand()}')

ELAPSED_TIME=$(($timestamp2 - $timestamp1))
echo "ELAPSED_TIME in second : $ELAPSED_TIME .Or in minute  :  $(($ELAPSED_TIME/60)) "
#kill -9 $client_pid



