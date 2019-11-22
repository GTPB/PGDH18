#!/bin/bash

myProjectName=$1
myInputFromSINS="../SINS_Sim/SINS_Output/$myProjectName/"
myOutputFolder="$myProjectName"
numberOfSimulations=$2
sinsSamplerDistributionJar="./SinsSampler.jar"

#dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
#move to current directory (where scrip is located)
cd ${0%/*}

java -jar $sinsSamplerDistributionJar -inputFromSINS $myInputFromSINS -numberOfSimulation $numberOfSimulations -hasClusters true -input $myProjectName -output $myOutputFolder -simulationName $myProjectName 

