import sys
import os
import argparse
import random
from math import *

with open("command.txt", "w") as of:
    of.write(" ".join(["python"]+sys.argv))

'''
This scripts runs hadd on single crystal files to 
group them in strips reading a DOF file
'''
parser = argparse.ArgumentParser()

#parser.add_argument("-f", "--files", type=str, help="input file", required=True)
parser.add_argument("-i", "--inputdir", type=str, help="Inputdir", required=True)
parser.add_argument("-o", "--outputdir", type=str, help="Outputdir", required=True)
parser.add_argument("-q", "--queue", type=str, help="Condor queue", default="longlunch", required=True)
parser.add_argument("--redo", action="store_true", default=False, help="Redo all files")
args = parser.parse_args()


# Prepare condor jobs
condor = '''executable              = run_script_mult.sh
output                  = output/strips.$(ClusterId).$(ProcId).out
error                   = error/strips.$(ClusterId).$(ProcId).err
log                     = log/strips.$(ClusterId).log
transfer_input_files    = run_script_mult.sh
+JobFlavour             = "{queue}"
+AccountingGroup	    = "group_u_CMS.CAF.ALCA"
queue arguments from arguments_mult.txt
'''

condor = condor.replace("{queue}", args.queue)
user = os.environ["USER"]

script = '''#!/bin/sh -e
export X509_USER_PROXY=/afs/cern.ch/user/{user1}/{user}/x509up_u35923
voms-proxy-info

JOBID=$1;  
INPUTFILE=$2;
OUTPUTFILE=$3;

xrdcp root://eoscms.cern.ch//eos/cms/store/user/lzygala/GJet/GJet_Validation_batch_mult.C .

root -l -b -q "GJet_Validation_batch.C++(\"${INPUTFILE}\",\"output_${JOBID}.root\")"

xrdcp ./output_${JOBID}.root root://eoscms.cern.ch/${OUTPUTFILE}

rm *.root
rm *.C
echo -e "DONE";
'''

script = script.replace("{user1}", user[:1])
script = script.replace("{user}", user)

arguments= []
if not os.path.exists(args.outputdir):
    os.makedirs(args.outputdir)

outputfiles = [args.outputdir +"/"+f for f in os.listdir(args.outputdir)]
inputfiles = [ f for f in os.listdir(args.inputdir)]

jobid = 0
for ifile in inputfiles:
    jobid +=1
    inputfile = args.inputdir + "/" + ifile
    outputfile = args.outputdir + "/" + ifile[:-5] + "_plots.root"

    if not args.redo and outputfile in outputfiles:
        continue
    
    arguments.append("{} {} {}".format(jobid,inputfile,outputfile))

print("Njobs: ", len(arguments))
    
with open("condor_job_mult.txt", "w") as cnd_out:
    cnd_out.write(condor)

with open("arguments_mult.txt", "w") as args:
    args.write("\n".join(arguments))

with open("run_script_multc.sh", "w") as rs:
    rs.write(script)

#os.system("condor_submit condor_job.txt")
