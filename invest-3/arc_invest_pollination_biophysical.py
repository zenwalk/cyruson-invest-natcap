#This script automatically generated by create_arc_scripts.py
import subprocess
import os
try:
    os.chdir('invest-3')
except:
    pass
subprocess.check_call('invest_pollination_biophysical.exe')
