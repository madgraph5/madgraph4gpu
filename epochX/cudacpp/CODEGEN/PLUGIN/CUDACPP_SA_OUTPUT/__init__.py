# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: O. Mattelaer (Sep 2021) for the MG5aMC CUDACPP plugin.
# Further modified by: O. Mattelaer, A. Valassi (2021-2023) for the MG5aMC CUDACPP plugin.

print('Load PLUGIN.CUDACPP_SA_OUTPUT')

# AV - Require Python >= 3.8 to ensure that {} dictionaries preserve the order of item insertion
# (note: python3.7 would probably be enough but this plugin has only been tested using python3.8)
import sys
minpython = (3,8)
if sys.version_info < minpython :

    print('ERROR! Cannot load PLUGIN.CUDACPP_SA_OUTPUT: Python >= %s.%s is required' % minpython)

else:

    # Import the required files
    # Example: import maddm_interface as maddm_interface # local file
    #          import madgraph.various.cluster as cluster # MG5 distribution file

    # Three types of functionality are allowed in a plugin
    #   1. new output mode
    #   2. new cluster support
    #   3. new interface

    # 1. Define new output mode.
    #    Example: new_output = {'myformat': MYCLASS}
    #    allows the command "output myformat PATH" in madgraph.
    #    MYCLASS should inherit from class madgraph.iolibs.export_v4.VirtualExporter
    import PLUGIN.CUDACPP_SA_OUTPUT.output as output
    new_output = { 'standalone_cudacpp' : output.PLUGIN_ProcessExporter }

    # 2. Define new way to handle the cluster.
    #    Example: new_cluster = {'mycluster': MYCLUSTERCLASS}
    #    allows the command "set cluster_type mycluster" in madgraph
    #    MYCLUSTERCLASS should inherit from class madgraph.various.cluster.Cluster.
    new_cluster = {}

    # 3. Define a new interface (allows adding/modifying MG5 command).
    #    This can be activated via ./bin/mg5_aMC --mode=PLUGINNAME.
    #    Put None if no dedicated command are required
    new_interface = None

    ########################## CONTROL VARIABLE ####################################

    __author__ = 'Andrea Valassi'
    __email__ = 'andrea.valassi@cern.ch'
    __version__ = (1,0,0)
    minimal_mg5amcnlo_version = (3,5,0)
    maximal_mg5amcnlo_version = (1000,1000,1000)
    latest_validated_version = (3,5,0)
