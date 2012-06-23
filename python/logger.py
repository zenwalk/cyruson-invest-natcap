from urllib import urlencode
from urllib2 import Request
from urllib2 import urlopen
import pkg_resources
import imp
import os.path

def log_model(model_name):
    """Submit a POST request to the defined URL with the modelname passed in as
    input.  The InVEST version number is also submitted, retrieved from the
    package's resources.

        model_name - a python string of the package version.

    returns nothing."""

    invest3_path = os.path.abspath(os.path.join(os.path.abspath(__file__),
        '..\\..\\invest-3\\invest_natcap\\__init__.pyc'))

    invest_natcap = imp.load_compiled('invest_natcap', invest3_path)
    invest_natcap.log_model(model_name + '_arc')
