# Copyright (C) 2020-2025 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: O. Mattelaer (Sep 2021) for the MG5aMC CUDACPP plugin.
# Further modified by: D. Massaro, O. Mattelaer, A. Valassi, Z. Wettersten (2021-2025) for the MG5aMC CUDACPP plugin.

# AV - Rename the plugin as CUDACPP_OUTPUT (even if the madgraph4gpu directory is still called CUDACPP_SA_OUTPUT)
# This can be used in mg5amcnlo in one of two ways:
# 1. production mode: a tarball containing CUDACPP_OUTPUT is untarred in PLUGIN/CUDACPP_OUTPUT
# 2. developers mode: the madgraph4gpu CUDACPP_SA_OUTPUT is symlinked in MG5aMC_PLUGIN/CUDACPP_OUTPUT (PYTHONPATH=.. ./bin/mg5_aMC -m CUDACPP_OUTPUT)
PLUGIN_NAME = __name__ # PLUGIN_NAME can be one of PLUGIN/CUDACPP_OUTPUT or MG5aMC_PLUGIN/CUDACPP_OUTPUT
print('Loading plugin %s'%PLUGIN_NAME)

# AV - Require Python >= 3.8 to ensure that {} dictionaries preserve the order of item insertion
# (note: python3.7 would probably be enough but this plugin has only been tested using python3.8)
import sys
minpython = (3,8)
if sys.version_info < minpython :

    print('ERROR! Cannot load plugin %s: Python >= %s.%s is required' % (PLUGIN_NAME, minpython[0], minpython[1] ))

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
    ###import PLUGIN.CUDACPP_OUTPUT.output as output # AV modify this to also allow MG5aMC_PLUGIN
    __import__('%s.output'%PLUGIN_NAME)
    output = sys.modules['%s.output'%PLUGIN_NAME]
    __import__('%s.trex'%PLUGIN_NAME)
    trex = sys.modules['%s.trex'%PLUGIN_NAME]
    new_output = { 'madevent_simd' : output.SIMD_ProcessExporter,
                   'madevent_gpu' : output.GPU_ProcessExporter,
                   'standalone_cudacpp' : output.PLUGIN_ProcessExporter,
                   'standalone_trex' : trex.TREX_ProcessExporter,
                   # the following one are used for the second exporter class 
                   # (not really needed so far but interesting if need
                   #  specialization in the futur) 
                   'standalone_simd' :  output.SIMD_ProcessExporter,
                   'standalone_cuda' :  output.GPU_ProcessExporter,
                  }
    new_reweight = {'madtrex': trex.TREX_ReweightInterface}

    # 2. Define new way to handle the cluster.
    #    Example: new_cluster = {'mycluster': MYCLUSTERCLASS}
    #    allows the command "set cluster_type mycluster" in madgraph
    #    MYCLUSTERCLASS should inherit from class madgraph.various.cluster.Cluster.
    new_cluster = {}

    # 3. Define a new interface (allows adding/modifying MG5 command).
    #    This can be activated via ./bin/mg5_aMC --mode=PLUGINNAME.
    #    Put None if no dedicated command are required
    if PLUGIN_NAME.rsplit('.',1)[0] == 'MG5aMC_PLUGIN':
        import madgraph.interface.master_interface as interface
        new_interface = interface.MasterCmd # use the default interface (but this is needed in the '-m' aka '--mode' option)
    else:
        new_interface = None

    ########################## CONTROL VARIABLE ####################################

    __author__ = 'Andrea Valassi'
    __email__ = 'andrea.valassi@cern.ch'

    # Plugin version (major,minor,patch) where major>1, 0<=minor<=99 and 0<=patch<=99
    # The release infrastructure expects 'vN.NN.NN' tags with 1-digit major and 2-digit minor and patch versions
    # and it takes care of converting the python tuple '(1,0,1)' into a version string 'v1.00.01'
    # NB! Do not use '(1,00,01)' here: leading zeros in decimal integer literals are not permitted in python (#1013)
    __version__ = (1,1,1)

    minimal_mg5amcnlo_version = (3,6,4)
    maximal_mg5amcnlo_version = (1000,1000,1000)
    latest_validated_version = (3,6,5)
