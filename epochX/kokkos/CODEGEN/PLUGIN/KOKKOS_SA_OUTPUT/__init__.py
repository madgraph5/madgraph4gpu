# AV - Require Python >= 3.8 to ensure that {} dictionaries preserve the order of item insertion
# (note: python3.7 would probably be enough but I test my scripts using python3.8)
import sys
minpython = (3, 8)
if sys.version_info < minpython: sys.exit("ERROR! Python >= %s.%s is required" % minpython)

## Import the required files
# Example: import maddm_interface as maddm_interface # local file
#          import madgraph.various.cluster as cluster #MG5 distribution file
# Three types of functionality are allowed in a plugin
#   1. new output mode
#   2. new cluster support
#   3. new interface

# 1. Define new output mode
#    example: new_output = {'myformat': MYCLASS}
#    madgraph will then allow the command "output myformat PATH"
#    MYCLASS should inherated of the class madgraph.iolibs.export_v4.VirtualExporter
import PLUGIN.KOKKOS_SA_OUTPUT.output as output
new_output = {'standalone_kokkos':output.PLUGIN_ProcessExporter}

# 2. Define new way to handle the cluster.
#    example new_cluster = {'mycluster': MYCLUSTERCLASS}
#    allow "set cluster_type mycluster" in madgraph
#    MYCLUSTERCLASS should inherated from madgraph.various.cluster.Cluster
new_cluster = {}

# 3. Define a new interface (allow to add/modify MG5 command)
#    This can be activated via ./bin/mg5_aMC --mode=PLUGINNAME
## Put None if no dedicated command are required
new_interface = None
 
########################## CONTROL VARIABLE ####################################

__author__ = ''
__email__ = ''
__version__ = (1,0,0)
minimal_mg5amcnlo_version = (2,9,5) 
maximal_mg5amcnlo_version = (1000,1000,1000)
latest_validated_version = (3,1,1)

