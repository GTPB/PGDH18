#!/bin/bash


myProjectName=$1
numberOfSimulations=1



nohup sh -c "./SINS_Sim/run_SINS.sh $myProjectName $numberOfSimulations >log_output/$myProjectName.SINS.log 2>&1 && ./SINS_Sampler/run_SINS_Sampler.sh $myProjectName $numberOfSimulations >log_output/$myProjectName.SAMPLER.log 2>&1" &


echo $! > ./save_pid/pid_$myProjectName.savePID
