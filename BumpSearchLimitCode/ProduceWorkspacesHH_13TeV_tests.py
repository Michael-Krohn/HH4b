
import os

masses =[]


for mass in masses:
  outputname = "submit_HH_Graviton_"+str(mass)+".src"
  logname = "submit_HH_Graviton_"+str(mass)+".log"
  outputfile = open(outputname,'w')
  outputfile.write('#!/bin/bash\n')
  outputfile.write("cd ${CMSSW_BASE}/src/HH4b/BumpSearchLimitCode_unblind_tau21cut_76X_May; eval `scramv1 run -sh`\n")

  #outputfile.write("root -b -q 'R2JJFitterHH_13TeV.cc("+str(mass)+","+'"Graviton"'+",true,50000.)'\n")


  outputfile.close()
  
  command="rm "+logname
  print command
  os.system(command)
  #command="bsub -q 1nh -o "+logname+" source "+outputname
  command="chmod 755 ./"+outputname+";./"+outputname
  print command
  os.system(command)

#masses=[1200, 1400, 1600, 1800, 2000, 2500]
masses=[2000]

for mass in masses:
  outputname = "submit_HH_Graviton_"+str(mass)+".src"
  logname = "submit_HH_Graviton_"+str(mass)+".log"
  outputfile = open(outputname,'w')
  outputfile.write('#!/bin/bash\n')
#  outputfile.write("cd ${CMSSW_BASE}/src/HH4b/BumpSearchLimitCode_unblind_tau21cut_76X_May; eval `scramv1 run -sh`\n")
  outputfile.write("eval `scramv1 run -sh`\n")
  outputfile.write("root -b -q 'R2JJFitterHH_13TeV.cc("+str(mass)+","+'"Graviton_subtr"'+",true,50000.)'\n")
  outputfile.write("root -b -q 'R2JJDatacardsMakerHH_13TeV.C("+str(mass)+","+'"Graviton_subtr"'+",true,50000.)'\n")
  outputfile.close()
  
  command="rm "+logname
  print command
  os.system(command)
  #command="bsub -q 1nh -o "+logname+" source "+outputname
  command="chmod 755 ./"+outputname+";./"+outputname
  print command
  os.system(command)
