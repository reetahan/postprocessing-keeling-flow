#!/usr/bin/env python
# coding: utf-8

# How to run this code? For scenario 0000 to 0099
# python get_nc_new.py 3586369 0 99

# ref: https://janakiev.com/blog/python-shell-commands/
import os
import sys
from tqdm import tqdm


# define the parameters for the command
job_id = sys.argv[1]
#csv_name = sys.argv[2]
sc_number_start = int(sys.argv[2])
sc_number_end = int(sys.argv[3])
# job_id = "10642359"
# #csv_name = "test"
# sc_number_start=0
# sc_number_end=2

# define the original directory
oridir = os.getcwd()
# define the filename 
filename = "run_new.sh"
# prefix of the destination
p1 = "/data/keeling/a/zzheng25/d/partmc-mam4-high-lat/"
p2 = "/scenarios/scenario_"


for i in tqdm(range(sc_number_start, sc_number_end+1)):
    sc_id = "%04i" %i  #define the file number
#     print("scenario id:", sc_id)
    
    dest = p1+job_id+p2+sc_id+"/"
    dest_file = dest+filename

    # define the command
    cpcmd = "cp"+" "+oridir+"/"+filename+" "+dest_file #copy
    # cdcmd = "cd" # change the directory
    execmd = "bash"+" "+filename # run the bash command
#     print(cpcmd)
#     print("\n")
#     print(dest)
#     print("\n")
#     print(execmd)
#     print("\n")

    # run the command
    os.system(cpcmd) # copy
    os.chdir(dest) # change the directory
    os.system(execmd)
#     print("Done")
