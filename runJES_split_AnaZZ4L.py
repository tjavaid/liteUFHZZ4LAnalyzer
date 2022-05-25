import os, time

if (not os.path.exists("Dump")):
    os.system("mkdir Dump")
if (not os.path.exists("Ntuples")):
    os.system("mkdir Ntuples")

dirMC = '/eos/cms/store/group/phys_muon/TagAndProbe/HZZ4L/2017/'
samplesMC  = [
'GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_2',
]



#njobs = 6
njobs = 1
for job in range(1,njobs+1):

  #if (not (job==2)): continue
  
  for sample in samplesMC:
    cmd = 'nohup ./ZZ4L_Ana.exe '+dirMC+'/'+sample+' Ntuples/'+sample+' 0 '+str(job)+' '+str(njobs)+' >& Dump/'+sample+'_'+str(job)+'.log &'
    print cmd
    os.system(cmd)


