# Adding logging to ncp-dev for this script. This is automatically generated by insert_logging.sh
import logger
logger.log_model('HRA_LaunchGUI')

import os, sys

if __name__=="__main__":
    os.system("start "+sys.argv[0][:sys.argv[0].rfind("\\")]+"\\HRA_SurveyTool.py")
    