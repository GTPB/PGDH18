#!/bin/bash


myProjectName=$1
numberOfSimulations=1



./SINS_Sim/run_SINS.sh $myProjectName $numberOfSimulations && ./SINS_Sampler/run_SINS_Sampler.sh $myProjectName $numberOfSimulations


