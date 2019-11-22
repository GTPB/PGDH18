nohup java -jar SINS2.jar -projectName 00MyTest -numberOfSimulation 10 -parallel yes -parCores 4 -outputFormat adegenet -outputFile ./output.txt -verbose false &

echo $! > save_pid 

tail -f nohup.out

