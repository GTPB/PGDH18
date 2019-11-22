#!/bin/bash


myProjectName=$1
myOutputFolder="./SINS_Output/"
numberOfSimulations=$2
sinsDistributionJar="./SINS2.jar"

#get current directory (where script is located)
#dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
#move to current directory (where script is located)
cd ${0%/*}

java -jar $sinsDistributionJar -projectName $myProjectName -compFormat noComp -numberOfSimulation $numberOfSimulations -takeSampledParametersFromFile no -outputFile $myOutputFolder -parallel no -v true
