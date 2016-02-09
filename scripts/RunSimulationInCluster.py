############################################################
# RunSimulationInCluster.py created by Javier Caravaca
#
# This script runs a RAT-PAC simulation divided in multiple
# different jobs that are run in parallel by the cluster
#
############################################################

import subprocess
import string
import os
import sys

def main():

    #Parse arguments
    if ( len(sys.argv) < 2):
        print 'Usage: python RunSimulationInCluster.py OUTPUTFILE NEVENTS'
        sys.exit()
    outfilename = sys.argv[1]
    nevents = int(sys.argv[2])

    #Max number of events per job
    neventsperjob = 100000
    njobs = nevents/neventsperjob

    print ' Number of events: ', nevents
    print ' Number of jobs: ', njobs

    #Create macro out of a template
    template_filename = '../mac/TheiaRnD_template.mac'
    template = string.Template(open(template_filename, 'r').read())

    for job in range(0,njobs):
        foutname = " \"../results/%s_%s.root\" " % (outfilename,job)
        rat_macro_text = template.substitute(OUTPUTFILE = foutname, NEVENTSPERJOB = neventsperjob)
        rat_macro_name = 'TheiaRnD_temp_%s.mac' % job
        log_name = 'TheiaRnD_temp_%s.log' % job
        rat_macro = file(rat_macro_name, 'w')
        rat_macro.write(rat_macro_text)
        rat_macro.close()
        #Run macro
        print 'Sending job %s to the cluster' % job
        subprocess.call("srun --comment=\" TheiaRnD\"  rat -q %s  > %s 2>&1 &" % (rat_macro_name, log_name) , shell=True)
        #Remove temporal macro file
        #os.remove(rat_macro_name)

if __name__ == "__main__":
    main()
