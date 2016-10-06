#!/bin/bash
file='mass.txt'
exec < $file 

while read line
do
    echo $line
    mv mass/"$line" MassPlotFineBins_subtr_Moriond_Silver.root
    cd ..
    root -l -q MiniTreeSignalProducerHH_13TeV.C++
    root -l -q MiniTreeProducerHH_13TeV.C++
    python ProduceWorkspacesHH_13TeV.py
    mv Graviton_subtr.root massOpt/"$line" 
    cd input
    mv MassPlotFineBins_subtr_Moriond_Silver.root mass/"$line"
done