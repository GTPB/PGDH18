#!/bin/bash


myProjectName=$1
sinsBuildInput="./SINS_BuildInput.jar"

#get current directory (where script is located)
#dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
#move to current directory (where script is located)
cd ${0%/*}

java -jar $sinsBuildInput -sinsIn ./SINS_Sim/input/$myProjectName -ssIn ./SINS_Sampler/input/$myProjectName

