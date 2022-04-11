import os
pjoin = os.path.join

import madgraph.iolibs.export_cpp as export_cpp

# AV - use template files from PLUGINDIR instead of MG5DIR
###from madgraph import MG5DIR
PLUGINDIR = os.path.dirname( __file__ )

# AV - model_handling includes custom UFOModelConverter and OneProcessExporter, plus additional patches
import PLUGIN.SYCL_SA_OUTPUT.model_handling as model_handling

#------------------------------------------------------------------------------------

# AV - modify misc.make_unique (remove a printout)
import madgraph.various.misc as misc

printordering = True
def PLUGIN_make_unique(input, keepordering=None):
    "remove duplicate in a list "
    global printordering
    if keepordering is None:
        keepordering = misc.madgraph.ordering
        if printordering:
            printordering = False
            misc.sprint('keepordering (default): %s'%keepordering) # AV - add a printout only in the first call
    else:
        misc.sprint('keepordering (argument): %s'%keepordering) # AV - add a printout at every call only if it is an argument	
    ###sprint(keepordering) # AV - remove the printout at every call
    if not keepordering:
        return list(set(input))
    else:
        return list(dict.fromkeys(input)) 

DEFAULT_make_unique = misc.make_unique
misc.make_unique = PLUGIN_make_unique

#------------------------------------------------------------------------------------

# AV - modify madgraph.iolibs.files.cp (preserve symlinks)
def PLUGIN_cp(path1, path2, log=True, error=False):
    """ simple cp taking linux or mix entry"""
    from madgraph.iolibs.files import format_path
    path1 = format_path(path1)
    path2 = format_path(path2)
    try:
        import shutil
        ###shutil.copy(path1, path2)
        shutil.copy(path1, path2, follow_symlinks=False) # AV
    except:
        from madgraph.iolibs.files import cp
        cp(path1, path2, log=log, error=error)

DEFAULT_cp = export_cpp.cp
export_cpp.cp = PLUGIN_cp

#------------------------------------------------------------------------------------

class PLUGIN_ProcessExporter(export_cpp.ProcessExporterGPU):
    # Class structure information
    #  - object
    #  - VirtualExporter(object) [in madgraph/iolibs/export_v4.py]
    #  - ProcessExporterCPP(VirtualExporter) [in madgraph/iolibs/export_cpp.py]
    #  - ProcessExporterGPU(ProcessExporterCPP) [in madgraph/iolibs/export_cpp.py]
    #      Note: only change class attribute
    #  - PLUGIN_ProcessExporter(ProcessExporterGPU)
    #      This class
    
    # Below are the class variable that are defined in export_v4.VirtualExporter
    # AV - keep defaults from export_v4.VirtualExporter
    # Check status of the directory. Remove it if already exists
    ###check = True 
    # Output type: [Template/dir/None] copy the Template (via copy_template), just create dir or do nothing 
    ###output = 'Template'

    # If sa_symmetry is true, generate fewer matrix elements
    # AV - keep OM's default for this plugin (using grouped_mode=False, "can decide to merge uu~ and u~u anyway")
    sa_symmetry = True

    # Below are the class variable that are defined in export_cpp.ProcessExporterGPU
    # AV - keep defaults from export_cpp.ProcessExporterGPU
    # Decide which type of merging is used [madevent/madweight]
    ###grouped_mode = False 
    # Other options
    ###default_opt = {'clean': False, 'complex_mass':False, 'export_format':'madevent', 'mp': False, 'v5_model': True }
    
    # AV - keep defaults from export_cpp.ProcessExporterGPU
    # AV - used in MadGraphCmd.do_output to assign export_cpp.ExportCPPFactory to MadGraphCmd._curr_exporter (if cpp or gpu)
    # AV - used in MadGraphCmd.export to assign helas_call_writers.(CPPUFO|GPUFO)HelasCallWriter to MadGraphCmd._curr_helas_model (if cpp or gpu)
    # Language type: 'v4' for f77, 'cpp' for C++ output
    exporter = 'gpu'

    # AV - use a custom OneProcessExporter
    ###oneprocessclass = export_cpp.OneProcessExporterGPU # responsible for P directory
    oneprocessclass = model_handling.PLUGIN_OneProcessExporter
    
    # Information to find the template file that we want to include from madgraph
    # you can include additional file from the plugin directory as well
    # AV - use template files from PLUGINDIR instead of MG5DIR and add gpu/mgOnGpuVectors.h
    # [NB: mgOnGpuConfig.h and check_sa.cu are handled through dedicated methods]
    ###s = MG5DIR + '/madgraph/iolibs/template_files/'
    s = PLUGINDIR + '/madgraph/iolibs/template_files/'
    from_template = {'src': [s+'gpu/rambo.h', s+'read_slha.h', s+'read_slha.cc',
                             s+'gpu/mgOnGpuTypes.h', s+'gpu/mgOnGpuVectors.h', s+'gpu/extras.h'],
                    'SubProcesses': [s+'gpu/timer.h', s+'gpu/timermap.h', s+'gpu/Memory.h', 
                                     s+'gpu/Makefile', s+'gpu/runTest.cc', s+'gpu/testxxx.cc', s+'gpu/testxxx_cc_ref.txt',
                                     s+'gpu/perf.py', s+'gpu/profile.sh']}
    to_link_in_P = ['timer.h', 'timermap.h', 'Memory.h', 'Makefile', 'runTest.cc', 'testxxx.cc', 'testxxx_cc_ref.txt', 'perf.py', 'profile.sh']

    # AV - use template files from PLUGINDIR instead of MG5DIR
    ###template_src_make = pjoin(MG5DIR, 'madgraph' ,'iolibs', 'template_files','gpu','Makefile_src')
    ###template_Sub_make = pjoin(MG5DIR, 'madgraph', 'iolibs', 'template_files','gpu','Makefile')
    template_src_make = pjoin(PLUGINDIR, 'madgraph' ,'iolibs', 'template_files','gpu','Makefile_src')
    template_Sub_make = pjoin(PLUGINDIR, 'madgraph', 'iolibs', 'template_files','gpu','Makefile')

    # AV - use a custom UFOModelConverter (model/aloha exporter)
    ###create_model_class =  export_cpp.UFOModelConverterGPU
    import PLUGIN.SYCL_SA_OUTPUT.model_handling as model_handling 
    create_model_class = model_handling.PLUGIN_UFOModelConverter
    
    # AV - "aloha_exporter" is not used anywhere!
    # (OM: "typically not defined but useful for this tutorial - the class for writing helas routine")
    ###aloha_exporter = None
    ###aloha_exporter = model_handling.PLUGIN_UFOHelasCallWriter

    # AV (default from OM's tutorial) - add a debug printout
    def __init__(self, *args, **kwargs):
        misc.sprint('Entering PLUGIN_ProcessExporter.__init__ (initialise the exporter)')
        return super().__init__(*args, **kwargs)

    # AV (default from OM's tutorial) - add a debug printout
    def copy_template(self, model):
        misc.sprint('Entering PLUGIN_ProcessExporter.copy_template (initialise the directory)')
        return super().copy_template(model)

    # AV - add debug printouts (in addition to the default one from OM's tutorial)
    def generate_subprocess_directory(self, subproc_group, fortran_model, me=None):
        misc.sprint('Entering PLUGIN_ProcessExporter.generate_subprocess_directory (create the directory)')
        misc.sprint('  type(subproc_group)=%s'%type(subproc_group)) # e.g. madgraph.core.helas_objects.HelasMatrixElement
        misc.sprint('  type(fortran_model)=%s'%type(fortran_model)) # e.g. madgraph.iolibs.helas_call_writers.GPUFOHelasCallWriter
        return super().generate_subprocess_directory(subproc_group, fortran_model, me)

    # AV (default from OM's tutorial) - add a debug printout
    def convert_model(self, model, wanted_lorentz=[], wanted_coupling=[]):
        misc.sprint('Entering PLUGIN_ProcessExporter.convert_model (create the model)')
        return super().convert_model(model, wanted_lorentz, wanted_coupling)

    # AV (default from OM's tutorial) - add a debug printout
    def finalize(self, matrix_element, cmdhistory, MG5options, outputflag):
        """Typically creating jpeg/HTML output/ compilation/...
           cmdhistory is the list of command used so far.
           MG5options are all the options of the main interface
           outputflags is a list of options provided when doing the output command"""
        misc.sprint('Entering PLUGIN_ProcessExporter.finalize')
        return super().finalize(matrix_element, cmdhistory, MG5options, outputflag)

    # AV (default from OM's tutorial) - overload settings and add a debug printout
    def modify_grouping(self, matrix_element):
        """allow to modify the grouping (if grouping is in place)
            return two value:
            - True/False if the matrix_element was modified
            - the new(or old) matrix element"""
        # Irrelevant here since group_mode=False so this function is never called
        misc.sprint('Entering PLUGIN_ProcessExporter.modify_grouping')
        return False, matrix_element

#------------------------------------------------------------------------------------
