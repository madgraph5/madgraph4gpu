# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: O. Mattelaer (Sep 2021) for the MG5aMC CUDACPP plugin.
# Further modified by: O. Mattelaer, A. Valassi, Z. Wettersten (2021-2023) for the MG5aMC CUDACPP plugin.

import os

# AV - load an independent 2nd copy of the export_cpp module (as PLUGIN_export_cpp) and use that within the plugin (workaround for #341)
# See https://stackoverflow.com/a/11285504
###import madgraph.iolibs.export_cpp as export_cpp # 1st copy
######import madgraph.iolibs.export_cpp as PLUGIN_export_cpp # this is not enough to define an independent 2nd copy: id(export_cpp)==id(PLUGIN_export_cpp)
import sys
import importlib.util
SPEC_EXPORTCPP = importlib.util.find_spec('madgraph.iolibs.export_cpp')
PLUGIN_export_cpp = importlib.util.module_from_spec(SPEC_EXPORTCPP)
SPEC_EXPORTCPP.loader.exec_module(PLUGIN_export_cpp)
sys.modules['PLUGIN.CUDACPP_SA_OUTPUT.PLUGIN_export_cpp'] = PLUGIN_export_cpp # allow 'import PLUGIN.CUDACPP_SA_OUTPUT.PLUGIN_export_cpp' in model_handling.py
del SPEC_EXPORTCPP
###print('id(export_cpp)=%s'%id(export_cpp))
###print('id(PLUGIN_export_cpp)=%s'%id(PLUGIN_export_cpp))

# AV - use template files from PLUGINDIR instead of MG5DIR
###from madgraph import MG5DIR
PLUGINDIR = os.path.dirname( __file__ )

# AV - model_handling includes the custom FileWriter, ALOHAWriter, UFOModelConverter, OneProcessExporter and HelasCallWriter, plus additional patches
import PLUGIN.CUDACPP_SA_OUTPUT.model_handling as model_handling

# AV - create a plugin-specific logger
import logging
logger = logging.getLogger('madgraph.PLUGIN.CUDACPP_SA_OUTPUT.output')

#------------------------------------------------------------------------------------

from os.path import join as pjoin
import madgraph.various.misc as misc

# AV - define the plugin's process exporter
# (NB: this is the plugin's main class, enabled in the new_output dictionary in __init__.py)
class PLUGIN_ProcessExporter(PLUGIN_export_cpp.ProcessExporterGPU):
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
    ###oneprocessclass = PLUGIN_export_cpp.OneProcessExporterGPU # responsible for P directory
    oneprocessclass = model_handling.PLUGIN_OneProcessExporter

    # Information to find the template file that we want to include from madgraph
    # you can include additional file from the plugin directory as well
    # AV - use template files from PLUGINDIR instead of MG5DIR and add gpu/mgOnGpuVectors.h
    # [NB: mgOnGpuConfig.h, check_sa.cc and fcheck_sa.f are handled through dedicated methods]
    ###s = MG5DIR + '/madgraph/iolibs/template_files/'
    s = PLUGINDIR + '/madgraph/iolibs/template_files/'
    from_template = {'.': [s+'.clang-format', s+'CMake/CMakeLists.txt',
                           s+'COPYRIGHT', s+'COPYING', s+'COPYING.LESSER' ],
                     'CMake': [s+'CMake/Compilers.txt', s+'CMake/Platforms.txt', s+'CMake/Macros.txt'],
                     'src': [s+'gpu/rambo.h', s+'read_slha.h', s+'read_slha.cc',
                             s+'gpu/mgOnGpuFptypes.h', s+'gpu/mgOnGpuCxtypes.h', s+'gpu/mgOnGpuVectors.h',
                             s+'CMake/src/CMakeLists.txt'],
                     'SubProcesses': [s+'gpu/nvtx.h', s+'gpu/timer.h', s+'gpu/timermap.h',
                                      s+'gpu/ompnumthreads.h', s+'gpu/CudaRuntime.h',
                                      s+'gpu/MemoryAccessHelpers.h', s+'gpu/MemoryAccessVectors.h',
                                      s+'gpu/MemoryAccessMatrixElements.h', s+'gpu/MemoryAccessMomenta.h',
                                      s+'gpu/MemoryAccessRandomNumbers.h', s+'gpu/MemoryAccessWeights.h',
                                      s+'gpu/MemoryAccessAmplitudes.h', s+'gpu/MemoryAccessWavefunctions.h',
                                      s+'gpu/MemoryAccessGs.h', s+'gpu/MemoryAccessCouplingsFixed.h',
                                      s+'gpu/MemoryAccessNumerators.h', s+'gpu/MemoryAccessDenominators.h',
                                      s+'gpu/EventStatistics.h', s+'gpu/CommonRandomNumbers.h',
                                      s+'gpu/CrossSectionKernels.cc', s+'gpu/CrossSectionKernels.h',
                                      s+'gpu/MatrixElementKernels.cc', s+'gpu/MatrixElementKernels.h',
                                      s+'gpu/RamboSamplingKernels.cc', s+'gpu/RamboSamplingKernels.h',
                                      s+'gpu/RandomNumberKernels.h',
                                      s+'gpu/CommonRandomNumberKernel.cc', s+'gpu/CurandRandomNumberKernel.cc',
                                      s+'gpu/Bridge.h', s+'gpu/BridgeKernels.cc', s+'gpu/BridgeKernels.h',
                                      s+'gpu/fbridge.cc', s+'gpu/fbridge.inc', s+'gpu/fsampler.cc', s+'gpu/fsampler.inc',
                                      s+'gpu/MadgraphTest.h', s+'gpu/runTest.cc',
                                      s+'gpu/testmisc.cc', s+'gpu/testxxx_cc_ref.txt',
                                      s+'gpu/perf.py', s+'gpu/profile.sh',
                                      s+'CMake/SubProcesses/CMakeLists.txt'],
                     'test': [s+'gpu/cudacpp_test.mk']}
    to_link_in_P = ['nvtx.h', 'timer.h', 'timermap.h',
                    'ompnumthreads.h', 'CudaRuntime.h',
                    'MemoryAccessHelpers.h', 'MemoryAccessVectors.h',
                    'MemoryAccessMatrixElements.h', 'MemoryAccessMomenta.h',
                    'MemoryAccessRandomNumbers.h', 'MemoryAccessWeights.h',
                    'MemoryAccessAmplitudes.h', 'MemoryAccessWavefunctions.h',
                    'MemoryAccessGs.h', 'MemoryAccessCouplingsFixed.h',
                    'MemoryAccessNumerators.h', 'MemoryAccessDenominators.h',
                    'EventStatistics.h', 'CommonRandomNumbers.h',
                    'CrossSectionKernels.cc', 'CrossSectionKernels.h',
                    'MatrixElementKernels.cc', 'MatrixElementKernels.h',
                    'RamboSamplingKernels.cc', 'RamboSamplingKernels.h',
                    'RandomNumberKernels.h', 'CommonRandomNumberKernel.cc', 'CurandRandomNumberKernel.cc',
                    'Bridge.h', 'BridgeKernels.cc', 'BridgeKernels.h',
                    'fbridge.cc', 'fbridge.inc', 'fsampler.cc', 'fsampler.inc',
                    'MadgraphTest.h', 'runTest.cc',
                    'testmisc.cc', 'testxxx_cc_ref.txt',
                    'cudacpp.mk', # this is generated from a template in Subprocesses but we still link it in P1
                    'testxxx.cc', # this is generated from a template in Subprocesses but we still link it in P1
                    'MemoryBuffers.h', # this is generated from a template in Subprocesses but we still link it in P1
                    'MemoryAccessCouplings.h', # this is generated from a template in Subprocesses but we still link it in P1
                    'perf.py', 'profile.sh']

    # AV - use template files from PLUGINDIR instead of MG5DIR and change their names
    ###template_src_make = pjoin(MG5DIR, 'madgraph' ,'iolibs', 'template_files','gpu','Makefile_src')
    ###template_Sub_make = pjoin(MG5DIR, 'madgraph', 'iolibs', 'template_files','gpu','Makefile')
    template_src_make = pjoin(PLUGINDIR, 'madgraph' ,'iolibs', 'template_files','gpu','cudacpp_src.mk')
    template_Sub_make = pjoin(PLUGINDIR, 'madgraph', 'iolibs', 'template_files','gpu','cudacpp.mk')
    template_tst_make = pjoin(PLUGINDIR, 'madgraph', 'iolibs', 'template_files','gpu','cudacpp_test.mk')

    # AV - use a custom UFOModelConverter (model/aloha exporter)
    ###create_model_class =  PLUGIN_export_cpp.UFOModelConverterGPU
    create_model_class = model_handling.PLUGIN_UFOModelConverter

    # AV - use a custom GPUFOHelasCallWriter
    # (NB: use "helas_exporter" - see class MadGraphCmd in madgraph_interface.py - not "aloha_exporter" that is never used!)
    ###helas_exporter = None
    helas_exporter = model_handling.PLUGIN_GPUFOHelasCallWriter # this is one of the main fixes for issue #341!

    # AV (default from OM's tutorial) - add a debug printout
    def __init__(self, *args, **kwargs):
        self.in_madevent_mode = False
        misc.sprint('Entering PLUGIN_ProcessExporter.__init__ (initialise the exporter)')
        return super().__init__(*args, **kwargs)

    # AV - overload the default version: create CMake directory, do not create lib directory
    def copy_template(self, model):
        misc.sprint('Entering PLUGIN_ProcessExporter.copy_template (initialise the directory)')
        try: os.mkdir(self.dir_path)
        except os.error as error: logger.warning(error.strerror + ' ' + self.dir_path)
        with misc.chdir(self.dir_path):
            logger.info('Creating subdirectories in directory %s' % self.dir_path)
            for d in ['src', 'Cards', 'SubProcesses', 'CMake', 'test', 'test/ref']: # AV - added CMake, test, test/ref; removed lib
                try: os.mkdir(d)
                except os.error as error: logger.warning(error.strerror + ' ' + os.path.join(self.dir_path,d))
            # Write param_card
            open(os.path.join('Cards','param_card.dat'), 'w').write(model.write_param_card())
            # Copy files in various subdirectories
            for key in self.from_template:
                for f in self.from_template[key]:
                    PLUGIN_export_cpp.cp(f, key) # NB this assumes directory key exists...
            # Copy src makefile
            if self.template_src_make:
                makefile_src = self.read_template_file(self.template_src_make) % {'model': self.get_model_name(model.get('name'))}
                open(os.path.join('src', 'cudacpp_src.mk'), 'w').write(makefile_src)
            # Copy SubProcesses makefile
            if self.template_Sub_make:
                makefile = self.read_template_file(self.template_Sub_make) % {'model': self.get_model_name(model.get('name'))}
                open(os.path.join('SubProcesses', 'cudacpp.mk'), 'w').write(makefile)
            # Copy test makefile
            if self.template_tst_make:
                makefile_test = self.read_template_file(self.template_tst_make) % {'model': self.get_model_name(model.get('name'))}
                open(os.path.join('test', 'cudacpp_test.mk'), 'w').write(makefile_test)

    # AV - add debug printouts (in addition to the default one from OM's tutorial)
    def generate_subprocess_directory(self, subproc_group, fortran_model, me=None):
        misc.sprint('Entering PLUGIN_ProcessExporter.generate_subprocess_directory (create the directory)')
        misc.sprint('  type(subproc_group)=%s'%type(subproc_group)) # e.g. madgraph.core.helas_objects.HelasMatrixElement
        misc.sprint('  type(fortran_model)=%s'%type(fortran_model)) # e.g. madgraph.iolibs.helas_call_writers.GPUFOHelasCallWriter
        misc.sprint('  type(me)=%s me=%s'%(type(me) if me is not None else None, me)) # e.g. int
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
        misc.sprint('Entering PLUGIN_ProcessExporter.finalize', self.in_madevent_mode)
        if self.in_madevent_mode:
            self.add_input_for_banner()
            path = os.path.realpath(os.curdir + os.sep + 'PLUGIN' + os.sep + 'CUDACPP_SA_OUTPUT')
            os.system(path + os.sep + 'patchMad.sh ' + self.dir_path + ' 32 PROD')
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


    def add_input_for_banner(self):
        # Note only for madevent mode
        new_parameters = ["{'name':'cudacpp_backend', 'value':'CPP', 'include':False, 'hidden':False}"
            ]

        finput = open(pjoin(self.dir_path, 'bin', 'internal', 'plugin_run_card'), 'w')
        
        for entry in new_parameters:
            finput.write(entry)

        
    
#------------------------------------------------------------------------------------
