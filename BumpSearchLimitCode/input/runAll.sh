#!/bin/bash
for ((k=1; k<=5; k=k+1 ))
do
    for ((j=1; j<=11; j=j+1 ))
    do
	width=$((k*5+15+90+j*5))
	if [ $width -gt 166 ]; then
	    continue;
	fi
	if [ $width -lt 116 ]; then
            continue;
        fi
	bmin=$((j*5+90))
	echo MassPlotFineBins_subtr_Moriond_Silver"$bmin"to"$width".root
	mv mass/MassPlotFineBins_subtr_Moriond_Silver"$bmin"to"$width".root MassPlotFineBins_subtr_Moriond_Silver.root
	cd ..
	root -l -q MiniTreeSignalProducerHH_13TeV.C++
	root -l -q MiniTreeProducerHH_13TeV.C++
	python ProduceWorkspacesHH_13TeV.py
	mv Graviton_subtr.root massOpt/"$bmin"to"$width".root 
	cd input
	mv MassPlotFineBins_subtr_Moriond_Silver.root mass/MassPlotFineBins_subtr_Moriond_Silver"$bmin"to"$width".root
    done
done


