#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Script to run one rembo simulation.
'''
import subprocess as subp 
import sys, os, argparse, shutil, glob, datetime, json
from runner.ext.namelist import Namelist

# Shortcut to namelist fuctionality
nml = Namelist()  # you could use your own module for reading namelist


def run_rembo():
    '''Main subroutine to run rembo.'''

    ### Manage command-line arguments ############################

    # Initialize argument parser
    parser = argparse.ArgumentParser()

    # Add options
    parser.add_argument('-e','--exe',type=str,default='rembo',
        help='''Define the executable file to use here. Shortcuts:
rembo = librembo/bin/test_rembo.x
''')
    parser.add_argument('-r','--run',action="store_true",
        help='Run the executable after preparing the job?')
    parser.add_argument('-s','--submit',action="store_true",
        help='Run the executable after preparing the job by submitting to the queue?')
    parser.add_argument('-q','--qos',type=str, default='short',
        help='Name of the qos the job should be submitted to (priority,short,medium,long)')
    parser.add_argument('-w','--wall', type=int, default=12,
        help='Maximum wall time to allow for job (only for jobs submitted to queue)')
    parser.add_argument('--email', type=str, default='None',
        help='Email address to send job notifications from cluster')
    parser.add_argument('--group', type=str, default='anthroia',
        help='Email address to send job notifications from cluster')
    parser.add_argument('-x',action="store_true",
        help='Use this argument if run_rembo.py is being called by job run')

    # Add arguments
    parser.add_argument('rundir',metavar='RUNDIR', type=str,
         help='Path where rembo simulation will run and store output.')
    parser.add_argument('par_path',metavar='PAR_PATH', type=str,
         help='Path to parameter file.')
    #parser.add_argument('par_path_2',metavar='PAR_PATH_2', type=str,
    #     help='Path to parameter file 2.',default='remboyelmo/yelmo_Greenland.nml')
    
    # Parse the arguments
    args = parser.parse_args()

    ### Manage user options and arguments ############################
    
    # Options
    exe_path    = args.exe       # Path relative to current working directory (cwd)
    run         = args.run 
    submit      = args.submit 
    qos         = args.qos  
    wtime       = args.wall 
    useremail   = args.email 
    usergroup   = args.group
    with_runner = args.x  

    # Arguments
    rundir      = args.rundir 
    par_path    = args.par_path    # Path relative to current working directory (cwd)
    #par_path_2  = args.par_path_2  # Path relative to current working directory (cwd)
    par_path_2  = "remboyelmo/yelmo_Greenland.nml"

    # Load simulation info from json configuration file 
    if os.path.isfile("run_config.json"):
        info = json.load(open("run_config.json"))
    else: 
        print("Required json file 'run_config.json' containing run options not found.")
        sys.exit()

    # Additional options, consistency checks

    # Copy the executable file to the output directory, or
    # call it from its compiled location?    
    copy_exec = True 

    # Submit overrides run 
    if submit: run = True 

    # Expand executable path shortcut if defined, otherwise exe_path remains unchanged.
    if exe_path in info["exe_shortcuts"]:
        exe_path = info["exe_shortcuts"].get(exe_path)

    # Also extract executable and path filenames 
    exe_fname   = os.path.basename(exe_path)
    par_fname   = os.path.basename(par_path)
    par_fname_2 = os.path.basename(par_path_2)

    # Check if using yelmo too
    if exe_fname == "yelmox_rembo.x":
        with_yelmo = True 
    else:
        with_yelmo = False 


    if with_yelmo:
        # Check yelmo related input files

        # Get path of constants parameter file
        const_path = "remboyelmo/yelmo_const_Earth.nml"

        # Make sure input files exist 
        if not os.path.isfile(const_path):
            print("Input file does not exist: {}".format(const_path))
            sys.exit() 

        if not os.path.isfile(par_path_2):
            print("Input file does not exist: {}".format(par_path_2))
            sys.exit() 

    if not os.path.isfile(par_path):
        print("Input file does not exist: {}".format(par_path))
        sys.exit() 

    if not os.path.isfile(exe_path):
        print("Input file does not exist: {}".format(exe_path))
        sys.exit() 
    
    ### Start the script to make the job, and then run it ############################


    # 1. Make the job (output directory, parameter files, additional data files/links)

    if with_runner:
        # Output directory already exists, update parameter file and write it to output directory
        runner_param_write(par_path,rundir)

    else: 
        # Make the output directory (and remove existing nc and nml files)
        # and copy the default parameter file 
        makedirs(rundir,remove=True)
        shutil.copy(par_path,rundir)
    
    # Add an output directory for 'obs' from model
    makedirs(rundir+"/obs",remove=True)

    if with_yelmo:
        # Copy the constants parameter file
        shutil.copy(const_path,rundir)
        shutil.copy(par_path_2,rundir)

    ## Generate symbolic links to input data folders
    for srcname in info["links"]:
        dstname = os.path.join(rundir,srcname)
        if os.path.islink(dstname): os.unlink(dstname)
        if os.path.islink(srcname):
            linkto = os.readlink(srcname)
            os.symlink(linkto, dstname)
        elif os.path.isdir(srcname):
            srcpath = os.path.abspath(srcname)
            os.symlink(srcpath,dstname)
        else:
            print("Warning: path does not exist {}".format(srcname))

    # Write the current git revision information to output directory 
    if os.path.isdir(".git"):
        head     = get_git_revision_hash()
        git_info = open(os.path.join(rundir,"git_revision"),'w').write(head)
    
    if with_yelmo:
        ## Generate symbolic link to maps folder
        srcname = "maps"
        dstname = os.path.join(rundir,srcname)
        srcpath = os.path.abspath(srcname)
        if os.path.islink(dstname): os.unlink(dstname)
        os.symlink(srcpath,dstname)

    # 2. Run the job

    # Generate the appropriate executable command to run job
    if copy_exec:
        # Assume executable is running from rundir
        executable = "./{}".format(exe_fname)

        # Also copy exe file to rundir 
        shutil.copy(exe_path,rundir)
        
    else:
        # Assume executable will run from current working directory 
        cwd = os.getcwd()
        executable = "{}/{}".format(cwd,exe_path)

    # Run the job if desired 
    if run:

        if submit:
            # Submit job to queue 
            pid = submitjob(rundir,executable,par_fname,qos,wtime,usergroup,useremail) 

        else:
            # Run job in background 
            pid = runjob(rundir,executable,par_fname)

    return 

######### Helper functions ############### 

def runjob(rundir,executable,par_path):
    '''Run a job generated with makejob.'''

    cmd = "cd {} && exec {} {} > {} &".format(rundir,executable,par_path,"out.out")

    print("Running job in background: {}".format(cmd))

    #os.system(cmd)
    proc = subp.Popen(cmd,shell=True,stdin=None,stdout=None,stderr=None,close_fds=True)
    #pid  = proc.pid+1   # This is not necessarily accurate - do not use for anything
    pid = 0

    # Alternative below is supposed to be more safe,
    # and provide the proper pid of the process itself,
    # but doesn't appear to actually work...
    # cmd = ['cd',rundir,'&&','exec',executable,par_path,'>','out.out &']
    # print " ".join(cmd)
    # proc = subp.Popen(cmd,shell=True,stdin=None,stdout=None,stderr=None,close_fds=True)
    # pid  = proc.pid
    #print "pid = {}".format(pid)

    return pid 

def submitjob(rundir,executable,par_path,qos,wtime,usergroup,useremail):
    '''Submit a job to a HPC queue (qsub,sbatch)'''

    # Get info about current system
    username  = os.environ.get('USER')
    hostname  = os.environ.get('HOSTNAME')

    # Command to be called 
    cmd = "{} {}".format(executable,par_path) 

    # Create the jobscript using current info
    nm_jobscript   = 'job.submit'
    path_jobscript = "{}/{}".format(rundir,nm_jobscript)
    
    if "cei" in hostname:
        script = jobscript_qsub(cmd,rundir,username,usergroup,wtime,useremail)
        jobfile = open(path_jobscript,'w').write(script)
        cmd_job = "cd {} && qsub {}".format(rundir,nm_jobscript)
        
    else:
        script  = jobscript_slurm(cmd,rundir,username,usergroup,qos,wtime,useremail)
        jobfile = open(path_jobscript,'w').write(script)
        cmd_job = "cd {} && sbatch {}".format(rundir,nm_jobscript)
    
    # Run the command (ie, change to output directory and submit job)
    #os.system(cmd_job)
    proc = subp.Popen(cmd_job,shell=True,stdin=None,stdout=None,stderr=None,close_fds=True)
    #pid  = proc.pid+1   # This is not necessarily accurate - do not use for anything
    pid = 0

    return pid 

def runner_param_write(par_path,rundir):
    '''Wrapper to perform parameter updates according to runner.json file 
       located in rundir, and then to write them to a new parameter file in rundir.
    '''

    # Read param file runner.json always written by runner to rundir
    rjson = os.path.join(rundir, 'runner.json')
    js    = json.load(open(rjson))

    # Load default namelist parameters
    params = nml.load(open(par_path))

    # Update default parameters with runner.json parameter values
    params.update(js['params']) 

    # To do: Add a map to be able to avoid defining the group name from the command line arguments
    #pmap      = {'p1':'g1.p1', 'p2':'g2.p2'}
    #nmlparams = {pmap.get(k,k):v for k,v in js['params'].items()}

    # To do: Maybe to automatically generate a map stripping the group name
    #pmap = {k.split('.')[-1]:k for k in params}

    # Write updated parameter file to rundir
    par_path_new = os.path.join(rundir,os.path.basename(par_path))
    nml.dump(params, open(par_path_new, 'w'))

    return 

def makedirs(dirname,remove):
    '''
    Make a directory (including sub-directories),
    but first ensuring that path doesn't already exist
    or some other error prevents the creation.
    '''

    try:
        os.makedirs(dirname)
        print('Directory created: {}'.format(dirname))
    except OSError:
        if os.path.isdir(dirname):
            print('Directory already exists: {}'.format(dirname))
            if remove:
                for f in glob.glob("{}/*.nc".format(dirname)): 
                    os.remove(f)
                for f in glob.glob("{}/*.nml".format(dirname)):
                    os.remove(f)
                #print('Removed *.nml and *.nc files.')
            pass
        else:
            # There was an error on creation, so make sure we know about it
            raise

    return

def get_git_revision_hash():
    #githash = subp.check_output(['git', 'describe', '--always', '--long', 'HEAD']).strip()
    githash = subp.check_output(['git', 'rev-parse', 'HEAD']).strip()
    return githash.decode("ascii") 

def jobscript_slurm(cmd,rundir,username,usergroup,qos,wtime,useremail):
    '''Definition of the job script'''

    jobname = "rembo" 
    
# Extra parameter options
##SBATCH --partition=ram_gpu
##SBATCH --mem=50000 
    
    # Check that wtime is consistent with qos
    if qos in ["priority","short"] and wtime > 24:
        print("Error in wtime for '{}'' queue, wtime = {}".format(qos,wtime))
        sys.exit()

    if qos == "medium" and wtime > 24*7:
        print("Error in wtime for '{}'' queue, wtime = {}".format(qos,wtime))
        sys.exit()
            
    script = """#! /bin/bash
#SBATCH --qos={}
#SBATCH --time={}:00:00
#SBATCH --job-name={}
#SBATCH --account={}
#SBATCH --mail-user={}
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --output=./out.out
#SBATCH --error=./out.err

# Run the job
{} 

""".format(qos,wtime,jobname,usergroup,useremail,cmd)

    return script

def jobscript_qsub(cmd,rundir,username,usergroup,wtime,useremail):
    '''Definition of the job script'''

    jobname = "rembo"
# Note: 
### jalv: this environment variable limits the number of CPUs used in eolo (=1, means one cpu)
###export OMP_NUM_THREADS=1
# Currently in .bashrc, but maybe should be done here.

    script = """#!/bin/bash
#$ -V                            # Ensure user enivronment variables are available
#$ -cwd                          # To use the current directory
#$ -m ae                         # Send mail when job is aborted (a), begins (b) and ends (e)
#$ -M {}                         # Send mail to this address
#$ -N  {}                        # (nombre del trabajo)
#$ -o ./out.out                  # (fichero de salida)   $JOB_ID
#$ -e ./out.err                  # (fichero de error)

#### Unused ####
####$ -l walltime={}:00:00       # Set wall time (hh:mm:ss)
####$ -l nodes=1:ppn=1           # Define as a serial job

# Run the job
time {}

""".format(useremail,jobname,wtime,cmd)

    return script


if __name__ == "__main__": 

    # Call main run_rembo function...
    run_rembo()



