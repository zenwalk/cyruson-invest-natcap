#This script automatically generated by create_arc_scripts.py
import sysimport osfrom distutils.sysconfig import get_python_libos.chdir(get_python_lib())import subprocess
os.chdir(sys.path[0])try:
    os.chdir('invest-3')
except:
    pass
subprocess.check_call('invest_pollination_valuation.exe')
