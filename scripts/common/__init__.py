import os
import sys

# enables to have "local" imports in scripts like primer_autodetection so that
# there is no need to copy and install the whole package to a docker file, but only the relevant scripts
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
