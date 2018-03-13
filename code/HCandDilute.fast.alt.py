#!/usr/bin/python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import datetime
import gzip
import os
import shutil
import subprocess
import sys
import time
import pybedtools

PD_DIST = 10
allelicPrimExe = "/mnt/webserver/datadisk/varcall/software/vcflib/bin/vcfallelicprimitives"

#-------------------------------------------------------------------------------------
# function to call subprocess.check_call
#-------------------------------------------------------------------------------------
def runAndTimeCmd(cmdToRun):
   timeStart = datetime.datetime.now()
   #print("Running cmd:" + cmdToRun)
   subprocess.check_call(cmdToRun, shell=True)
   timeEnd = datetime.datetime.now()
   #print("Time taken: " + str(timeEnd-timeStart))
   
#------------------------------------------------------------------------------------
# filter function to filter features present in hash. customized for mutect output. 
#-----------------------------------------------------------------------------------
def presentFilter_1(feature, vcfDict):
   if (feature[0], feature[1], feature[3]) in vcfDict:
      return False
   else:
      return True

#------------------------------------------------------------------------------------
# filter function to filter features present in hash. customized for variant calling output. 
#-----------------------------------------------------------------------------------
def presentFilter_2(feature, vcfDict):
   if (feature[0], feature[1], feature[2]) in vcfDict:
      return False
   else:
      return True

#------------------------------------------------------------------------------------
# subtract vcfB calls from vcfA calls: subtractBed is not able to get this right, customized for mutect output.
#------------------------------------------------------------------------------------
def subtractCalls_1(vcfA, vcfB, fileName):
   bCalls = {}
   # load calls from vcfB
   for feature in vcfB:
      bCalls[(feature[0], feature[1], feature[3])] = feature[4]
   outVcf = vcfA.filter(presentFilter_1, vcfDict=bCalls)
   outVcf.saveas(fileName)

#------------------------------------------------------------------------------------
# subtract vcfB calls from vcfA calls: subtractBed is not able to get this right, customized for variant calling output.
#------------------------------------------------------------------------------------
def subtractCalls_2(vcfA, vcfB, fileName):
   bCalls = {}
   # load calls from vcfB
   for feature in vcfB:
      bCalls[(feature[0], feature[1], feature[3])] = feature[4]
   outVcf = vcfA.filter(presentFilter_2, vcfDict=bCalls)
   outVcf.saveas(fileName)

#-------------------------------------------------------------------------------------
# limit candidate variants to panel and nist bed region
#-------------------------------------------------------------------------------------
def cleanVcf_2(runPath, vcfName, RunOutVcf, sharedRegion, overlap, dilutionVcfFile = None):
   os.chdir(runPath)
   sharedBed = pybedtools.BedTool(sharedRegion)

   # create vcf with no chr
   runAndTimeCmd("tail -n+2 " + vcfName + " | sed 's/chr//' > tmp.vcf")
   # add header
   runAndTimeCmd("echo '##fileformat=VCFv4.1' > header.txt") 
   runAndTimeCmd("cat header.txt tmp.vcf > nochr.vcf")

   # load vcfs and limit them to the shared regions
   runVcf = pybedtools.BedTool("nochr.vcf").intersect(sharedBed, f=overlap, u=True)

   # if dilution vcf is provided, subtract those mutations out. Note: dilution vcf must have been intersected with the shared region BED. 
   if dilutionVcfFile != None:
      dilutionVcf =  pybedtools.BedTool(dilutionVcfFile)
      dilutionCalls = len(dilutionVcf)
      subtractCalls_2(runVcf, dilutionVcf, RunOutVcf)
   else:
      runVcf.saveas(RunOutVcf)

   os.remove('nochr.vcf')
   os.remove('tmp.vcf')
   os.remove('header.txt')

#-------------------------------------------------------------------------------------
# main program for running from shell
#-------------------------------------------------------------------------------------
if __name__ == "__main__":
   runPath = sys.argv[1]
   vcfName = sys.argv[2]
   RunOutVcf = sys.argv[3]
   sharedRegion = sys.argv[4]
   overlap = float(sys.argv[5])

   if len(sys.argv) > 6 and sys.argv[6] != "None":
      dilutionVcfFile = sys.argv[6]
   else:
      dilutionVcfFile = None
   cleanVcf_2(runPath, vcfName, RunOutVcf, sharedRegion, overlap, dilutionVcfFile)
