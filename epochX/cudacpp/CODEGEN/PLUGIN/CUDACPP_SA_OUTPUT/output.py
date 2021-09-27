import madgraph.iolibs.export_cpp as export_cpp
import madgraph.various.misc as misc

from madgraph import MG5DIR

import os
pjoin = os.path.join

class MY_CPP_Standalone(export_cpp.ProcessExporterGPU):
    # class structure information
    # object
    #  - VirtualExporter(object) [in madgraph/iolibs/export_v4.py]
    #  - ProcessExporterCPP(VirtualExporter) [in madgraph/iolibs/export_cpp.py]
    #  - ProcessExporterGPU(ProcessExporterCPP) [in madgraph/iolibs/export_cpp.py]
    #      Note: only change class attribute
    #  - MY_CPP_Standalone(ProcessExporterGPU)
    #      This class

    
    # check status of the directory. Remove it if already exists
    check = True 
    # Language type: 'v4' for f77/ 'cpp' for C++ output
    exporter = 'gpu' 
    # Output type:
    #[Template/dir/None] copy the Template, just create dir  or do nothing 
    output = 'Template'
    # Decide which type of merging is used [madevent/madweight]
    grouped_mode = False
    # if no grouping on can decide to merge uu~ and u~u anyway:
    sa_symmetry = True


    # Here below all the class variable that are define for the GPU mode.
    # technically they are no need to overwrite those
    grouped_mode = False
    #exporter = 'gpu'
    
    default_opt = {'clean': False, 'complex_mass':False,
                        'export_format':'madevent', 'mp': False,
                        'v5_model': True
                        }
    
    ###oneprocessclass = export_cpp.OneProcessExporterGPU # responsible for P directory
    import PLUGIN.CUDACPP_SA_OUTPUT.model_handling as model_handling
    oneprocessclass = model_handling.OneProcessExporterGPU
    
    # information to find the template file that we want to include from madgraph
    # you can include additional file from the plugin directory as well
    PLUGINDIR = os.path.dirname( __file__ )
    ###s = MG5DIR + '/madgraph/iolibs/template_files/'
    s = PLUGINDIR + '/madgraph/iolibs/template_files/'
    from_template = {'src': [s+'gpu/rambo.h', s+'gpu/rambo.cc', s+'read_slha.h', s+'read_slha.cc',
                             s+'gpu/mgOnGpuTypes.h', s+'gpu/grambo.cu'],
                    'SubProcesses': [s+'gpu/timer.h', s+'gpu/Makefile', s+'gpu/nvtx.h',
                                     s+'gpu/nvtx.h', s+'gpu/check.cc',
                                     s+'gpu/timermap.h', s+'gpu/profile.sh',
                                     s+'gpu/perf.py', s+ 'gpu/Memory.h', s + 'gpu/runTest.cc']}
    to_link_in_P = ['Makefile', 'timer.h', 'timermap.h', 'nvtx.h', 'perf.py', 'Memory.h', 'runTest.cc']

    
    ###template_src_make = pjoin(MG5DIR, 'madgraph' ,'iolibs', 'template_files','gpu','Makefile_src')
    ###template_Sub_make = pjoin(MG5DIR, 'madgraph', 'iolibs', 'template_files','gpu','Makefile')
    template_src_make = pjoin(PLUGINDIR, 'madgraph' ,'iolibs', 'template_files','gpu','Makefile_src')
    template_Sub_make = pjoin(PLUGINDIR, 'madgraph', 'iolibs', 'template_files','gpu','Makefile')

    #For model/aloha exporter (typically not used)
    create_model_class =  export_cpp.UFOModelConverterGPU
    #import PLUGIN.CUDACPP_SA_OUTPUT.model_handling as model_handling 
    #create_model_class = model_handling.UFOModelConverterGPU
    
    # typically not defined but usufull for this tutorial the class for writing helas routine
    #aloha_exporter = None
    #aloha_exporter = model_handling.GPUFOHelasCallWriter
    

    def __init__(self, *args, **opts):
        misc.sprint("Initialise the exporter")
        return super(MY_CPP_Standalone, self).__init__(*args, **opts)

    def copy_template(self, model):

        misc.sprint("initialise the directory")
        return super(MY_CPP_Standalone, self).copy_template(model)


    def generate_subprocess_directory(self, subproc_group,
                                         fortran_model, me=None):
        
        misc.sprint('create the directory')
        return super(MY_CPP_Standalone, self).generate_subprocess_directory(subproc_group, fortran_model, me)


    def convert_model(self, model, wanted_lorentz=[], wanted_coupling=[]):
        misc.sprint('create the model')
        return super(MY_CPP_Standalone, self).convert_model(model, wanted_lorentz, wanted_coupling)

    def finalize(self, matrix_element, cmdhistory, MG5options, outputflag):
        """typically creating jpeg/HTML output/ compilation/...
           cmdhistory is the list of command used so far.
           MG5options are all the options of the main interface
           outputflags is a list of options provided when doing the output command"""
        misc.sprint("pass here")
        return super(MY_CPP_Standalone, self).finalize(matrix_element, cmdhistory, MG5options, outputflag)

    def modify_grouping(self, matrix_element):
        """allow to modify the grouping (if grouping is in place)
            return two value:
            - True/False if the matrix_element was modified
            - the new(or old) matrix element"""
        #irrelevant here since group_mode=False so this function is never called
        return False, matrix_element


    
