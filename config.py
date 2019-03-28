'''
Script to generate Makefile with the correct compiler
configuration for a given machine.
'''

from subprocess import *
import argparse
import sys 

# Manage command-line arguments to load compiler option
parser = argparse.ArgumentParser()
parser.add_argument(
    'compiler', metavar='COMPILER', type=str, choices=['ifort','gfortran'],
     help='Compiler choice should be included as an argument: ifort or gfortran')
args = parser.parse_args()


# Determine locations
config_dir = "config/"
target_dir = "./"
compiler   = args.compiler 

# Determine the current hostname
proc = Popen("hostname",shell=True,stdout=PIPE,stderr=PIPE)
host = proc.communicate()[0].strip().decode("utf-8")

if host == "manager.cei.ucm.local": host = "eolo"
if host in ["login01","login02"]:   host = "pik"

print("{}{}_{}".format(config_dir,host,compiler))

# Load template Makefile and insert compiler configuration
makefile0    = open(config_dir+"Makefile").read()
compile_info = open(config_dir+host+"_"+compiler).read()
makefile1 = makefile0.replace("<COMPILER_CONFIGURATION>",compile_info)

# Write the new Makefile to the main directory
open(target_dir+"Makefile","w").write(makefile1)

print( "".join(["\n\nMakefile configuration complete for host - compiler: ",host," - ",compiler,"\n"]) )

instructions = '''==== How to run rembo ====\n
# Make a link to the ice_data path, if available.
ln -s path/to/ice_data ice_data 

# Compile a test program
make clean 
make rembo      # makes executable librembo/bin/rembo.x 

# Run the program using run_rembo.py (see run_rembo.py -h for details) 
python run_rembo.py -r -e rembo output/test options_rembo

# Check the output 
cd output/test 
ncview clima.nc 
'''

print(instructions)


