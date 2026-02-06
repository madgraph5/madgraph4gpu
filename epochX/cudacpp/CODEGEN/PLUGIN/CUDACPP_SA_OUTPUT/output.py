# Copyright (C) 2020-2025 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Sep 2021) for the MG5aMC CUDACPP plugin.
# Further modified by: S. Hageboeck, O. Mattelaer, S. Roiser, J. Teig, A. Valassi, Z. Wettersten (2021-2024) for the MG5aMC CUDACPP plugin.

import os
import sys
import subprocess

# AV - PLUGIN_NAME can be one of PLUGIN/CUDACPP_OUTPUT or MG5aMC_PLUGIN/CUDACPP_OUTPUT
PLUGIN_NAME = __name__.rsplit('.',1)[0]

# AV - load an independent 2nd copy of the export_cpp module (as PLUGIN_export_cpp) and use that within the plugin (workaround for #341)
# See https://stackoverflow.com/a/11285504
###import madgraph.iolibs.export_cpp as export_cpp # 1st copy
######import madgraph.iolibs.export_cpp as PLUGIN_export_cpp # this is not enough to define an independent 2nd copy: id(export_cpp)==id(PLUGIN_export_cpp)
import importlib.util
SPEC_EXPORTCPP = importlib.util.find_spec('madgraph.iolibs.export_cpp')
PLUGIN_export_cpp = importlib.util.module_from_spec(SPEC_EXPORTCPP)
SPEC_EXPORTCPP.loader.exec_module(PLUGIN_export_cpp)
###sys.modules['PLUGIN.CUDACPP_OUTPUT.PLUGIN_export_cpp'] = PLUGIN_export_cpp # allow 'import PLUGIN.CUDACPP_OUTPUT.PLUGIN_export_cpp' in model_handling.py
sys.modules['%s.PLUGIN_export_cpp'%PLUGIN_NAME] = PLUGIN_export_cpp # allow 'import <PLUGIN_NAME>.PLUGIN_export_cpp' in model_handling.py
del SPEC_EXPORTCPP
###print('id(export_cpp)=%s'%id(export_cpp))
###print('id(PLUGIN_export_cpp)=%s'%id(PLUGIN_export_cpp))

# AV - use template files from PLUGINDIR instead of MG5DIR
###from madgraph import MG5DIR
PLUGINDIR = os.path.dirname( __file__ )

# AV - model_handling includes the custom FileWriter, ALOHAWriter, UFOModelConverter, OneProcessExporter and HelasCallWriter, plus additional patches
###import PLUGIN.CUDACPP_OUTPUT.model_handling as model_handling # AV modify this to also allow MG5aMC_PLUGIN
__import__('%s.model_handling'%PLUGIN_NAME)
model_handling = sys.modules['%s.model_handling'%PLUGIN_NAME]

# AV - create a plugin-specific logger
import logging
logger = logging.getLogger('madgraph.%s.output'%PLUGIN_NAME)
from madgraph import MG5DIR
#------------------------------------------------------------------------------------

from os.path import join as pjoin
import madgraph.iolibs.files as files
import madgraph.iolibs.export_v4 as export_v4
import madgraph.various.misc as misc

from . import launch_plugin


# AV - define the plugin's process exporter
# (NB: this is the plugin's main class, enabled in the new_output dictionary in __init__.py)
class PLUGIN_ProcessExporter(PLUGIN_export_cpp.ProcessExporterCPP):
    # Class structure information
    #  - object
    #  - VirtualExporter(object) [in madgraph/iolibs/export_v4.py]
    #  - ProcessExporterCPP(VirtualExporter) [in madgraph/iolibs/export_cpp.py]
    #  - PLUGIN_ProcessExporter(ProcessExporterCPP)
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
    grouped_mode = False
    # Other options
    default_opt = {'clean': False, 'complex_mass':False, 'export_format':'madevent', 'mp': False, 'v5_model': True }

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
                             s+'gpu/constexpr_math.h',
                             s+'gpu/cudacpp_config.mk',
                             s+'CMake/src/CMakeLists.txt' ],
                     'SubProcesses': [s+'gpu/nvtx.h', s+'gpu/timer.h', s+'gpu/timermap.h',
                                      s+'gpu/ompnumthreads.h', s+'gpu/GpuRuntime.h', s+'gpu/GpuAbstraction.h',
                                      s+'gpu/color_sum.h',
                                      s+'gpu/MemoryAccessHelpers.h', s+'gpu/MemoryAccessVectors.h',
                                      s+'gpu/MemoryAccessMatrixElements.h', s+'gpu/MemoryAccessMomenta.h',
                                      s+'gpu/MemoryAccessRandomNumbers.h', s+'gpu/MemoryAccessWeights.h',
                                      s+'gpu/MemoryAccessAmplitudes.h', s+'gpu/MemoryAccessWavefunctions.h',
                                      s+'gpu/MemoryAccessGs.h', s+'gpu/MemoryAccessCouplingsFixed.h',
                                      s+'gpu/MemoryAccessNumerators.h', s+'gpu/MemoryAccessDenominators.h',
                                      s+'gpu/MemoryAccessChannelIds.h', s+'gpu/MemoryAccessIflavorVec.h',
                                      s+'gpu/EventStatistics.h', s+'gpu/CommonRandomNumbers.h',
                                      s+'gpu/CrossSectionKernels.cc', s+'gpu/CrossSectionKernels.h',
                                      s+'gpu/MatrixElementKernels.cc', s+'gpu/MatrixElementKernels.h',
                                      s+'gpu/RamboSamplingKernels.cc', s+'gpu/RamboSamplingKernels.h',
                                      s+'gpu/RandomNumberKernels.h', s+'gpu/CommonRandomNumberKernel.cc',
                                      s+'gpu/CurandRandomNumberKernel.cc', s+'gpu/HiprandRandomNumberKernel.cc',
                                      s+'gpu/Bridge.h', s+'gpu/BridgeKernels.cc', s+'gpu/BridgeKernels.h',
                                      s+'gpu/fbridge.cc', s+'gpu/fbridge.h', s+'gpu/fbridge.inc', s+'gpu/fsampler.cc', s+'gpu/fsampler.inc',
                                      s+'gpu/MadgraphTest.h', s+'gpu/runTest.cc',
                                      s+'gpu/testmisc.cc', s+'gpu/testxxx_cc_ref.txt', s+'gpu/valgrind.h',
                                      s+'gpu/perf.py', s+'gpu/profile.sh',
                                      s+'gpu/cudacpp_overlay.mk', s+'gpu/makefile_wrapper.mk',
                                      s+'gpu/umami.h', s+'gpu/umami.cc',
                                      s+'CMake/SubProcesses/CMakeLists.txt'],
                     'test': [s+'gpu/cudacpp_test.mk']}

    to_link_in_P = ['nvtx.h', 'timer.h', 'timermap.h',
                    'ompnumthreads.h', 'GpuRuntime.h', 'GpuAbstraction.h',
                    'color_sum.h',
                    'MemoryAccessHelpers.h', 'MemoryAccessVectors.h',
                    'MemoryAccessMatrixElements.h', 'MemoryAccessMomenta.h',
                    'MemoryAccessRandomNumbers.h', 'MemoryAccessWeights.h',
                    'MemoryAccessAmplitudes.h', 'MemoryAccessWavefunctions.h',
                    'MemoryAccessGs.h', 'MemoryAccessCouplingsFixed.h',
                    'MemoryAccessNumerators.h', 'MemoryAccessDenominators.h',
                    'MemoryAccessChannelIds.h', 'MemoryAccessIflavorVec.h',
                    'EventStatistics.h', 'CommonRandomNumbers.h',
                    'CrossSectionKernels.cc', 'CrossSectionKernels.h',
                    'MatrixElementKernels.cc', 'MatrixElementKernels.h',
                    'RamboSamplingKernels.cc', 'RamboSamplingKernels.h',
                    'RandomNumberKernels.h', 'CommonRandomNumberKernel.cc',
                    'CurandRandomNumberKernel.cc', 'HiprandRandomNumberKernel.cc',
                    'Bridge.h', 'BridgeKernels.cc', 'BridgeKernels.h',
                    'fbridge.cc', 'fbridge.h', 'fbridge.inc', 'fsampler.cc', 'fsampler.inc',
                    'MadgraphTest.h', 'runTest.cc',
                    'testmisc.cc', 'testxxx_cc_ref.txt', 'valgrind.h',
                    'cudacpp.mk', # this is generated from a template in Subprocesses but we still link it in P1
                    'cudacpp_overlay.mk', # this is generated from a template in Subprocesses but we still link it in P1
                    'testxxx.cc', # this is generated from a template in Subprocesses but we still link it in P1
                    'MemoryBuffers.h', # this is generated from a template in Subprocesses but we still link it in P1
                    'MemoryAccessCouplings.h', # this is generated from a template in Subprocesses but we still link it in P1
                    'umami.h', 'umami.cc',
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
        self.in_madevent_mode = False # see MR #747
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

    # OM - overload export_v4.py version to add additional_clean section (and avoid patchMad.sh for Source/makefile)
    def write_source_makefile(self, writer, model=None, default=None):
        if default:
            replace_dict = default
        else:
            raise Exception('primary exporter should have been run first')
        path = pjoin(PLUGINDIR , 'madgraph', 'iolibs', 'template_files', 'madevent_makefile_source_addon')
        replace_dict['additional_clean'] += open(path).read()
        if writer:
            path = pjoin(MG5DIR, 'madgraph', 'iolibs','template_files','madevent_makefile_source')
            text = open(path).read() % replace_dict
            writer.write(text)

    # AV - add debug printouts (in addition to the default one from OM's tutorial)
    def generate_subprocess_directory(self, subproc_group, fortran_model, me=None):
        misc.sprint('Entering PLUGIN_ProcessExporter.generate_subprocess_directory (create the directory)')
        misc.sprint('  type(subproc_group)=%s'%type(subproc_group)) # e.g. madgraph.core.helas_objects.HelasMatrixElement
        misc.sprint('  type(fortran_model)=%s'%type(fortran_model)) # e.g. madgraph.iolibs.helas_call_writers.GPUFOHelasCallWriter
        misc.sprint('  type(me)=%s me=%s'%(type(me) if me is not None else None, me)) # e.g. int
        misc.sprint("need to link", self.to_link_in_P)
        out = super().generate_subprocess_directory(subproc_group, fortran_model, me)
        return out
    # AV (default from OM's tutorial) - add a debug printout
    def convert_model(self, model, wanted_lorentz=[], wanted_couplings=[]):
        if hasattr(model , 'cudacpp_wanted_ordered_couplings'):
            wanted_couplings = model.cudacpp_wanted_ordered_couplings
            del model.cudacpp_wanted_ordered_couplings
        return super().convert_model(model, wanted_lorentz, wanted_couplings)

    # AV (default from OM's tutorial) - add a debug printout
    def finalize(self, matrix_element, cmdhistory, MG5options, outputflag):
        """Typically creating jpeg/HTML output/ compilation/...
            cmdhistory is the list of command used so far.
            MG5options are all the options of the main interface
            outputflags is a list of options provided when doing the output command"""
        ###misc.sprint('Entering PLUGIN_ProcessExporter.finalize', self.in_madevent_mode, type(self))
        if self.in_madevent_mode:
            # Modify makefiles and symlinks to avoid doing
            # make -f makefile -f cudacpp_overlay.mk to include the overlay
            # and instead just use `make`, see #1052
            subprocesses_dir = pjoin(self.dir_path, "SubProcesses")
            files.cp(pjoin(subprocesses_dir, "makefile"), pjoin(subprocesses_dir, "makefile_original.mk"))
            files.rm(pjoin(subprocesses_dir, "makefile"))
            files.ln(pjoin(subprocesses_dir, "makefile_wrapper.mk"), subprocesses_dir, 'makefile')

            patch_coupl_write = r"""set -euo pipefail
# Get last fields from lines starting with WRITE(*,2)
gcs=$(awk '$1=="WRITE(*,2)" {print $NF}' coupl_write.inc)

for gc in $gcs; do
  if grep -q "$gc(VECSIZE_MEMMAX)" coupl.inc; then
    awk -v gc="$gc" '{
      if ($1=="WRITE(*,2)" && $NF==gc) print $0"(1)";
      else print
    }' coupl_write.inc > coupl_write.inc.new
    mv coupl_write.inc.new coupl_write.inc
  fi
done"""
            try:
                result = subprocess.run(
                    ["bash", "-c", patch_coupl_write],
                    cwd=pjoin(self.dir_path, "Source", "MODEL"),
                    text=True,
                    capture_output=True,
                    check=True,  # raise CalledProcessError on non-zero exit
                )
                misc.sprint(result.returncode)
            except subprocess.CalledProcessError as e:
                logger.debug("####### \n stdout is \n %s", e.stdout)
                logger.info("####### \n stderr is \n %s", e.stderr)
                logger.info("return code is %s\n", e.returncode)
                raise Exception("ERROR while patching coupl_write.inc") from e

            # Additional patching (OM)
            self.add_madevent_plugin_fct() # Added by OM
        # do not call standard finalize since is this is already done...
        #return super().finalize(matrix_element, cmdhistory, MG5options, outputflag)

    # AV (default from OM's tutorial) - overload settings and add a debug printout
    def modify_grouping(self, matrix_element):
        """allow to modify the grouping (if grouping is in place)
            return two value:
            - True/False if the matrix_element was modified
            - the new(or old) matrix element"""
        # Irrelevant here since group_mode=False so this function is never called
        misc.sprint('Entering PLUGIN_ProcessExporter.modify_grouping')
        return False, matrix_element

    # OM adding a new way to "patch" python file such that the launch command of MG5aMC is working
    # this consist in a file plugin_interface.py
    # which contains a series of functions and one dictionary variable TO_OVERWRITE
    # that will be used to have temporary overwrite of all the key variable passed as string by their value.
    # all variable that are file related should be called as madgraph.dir.file.variable
    def add_madevent_plugin_fct(self):
        """this consist in a file plugin_interface.py
        which contains a series of functions and one dictionary variable TO_OVERWRITE
        that will be used to have temporary overwrite of all the key variable passed as string by their value.
        all variable that are file related should be called as madgraph.dir.file.variable
        """
        plugin_path = os.path.dirname(os.path.realpath( __file__ ))
        files.cp(pjoin(plugin_path, 'launch_plugin.py'), pjoin(self.dir_path, 'bin', 'internal'))
        files.ln(pjoin(self.dir_path, 'lib'),  pjoin(self.dir_path, 'SubProcesses'))

#------------------------------------------------------------------------------------

class PLUGIN_ProcessExporter_MadEvent(PLUGIN_ProcessExporter):
    """ a class to include all tweak related to madevent and not related to standalone.
        in practise this class is never called but only the SIMD or GPU related class"""

    s = PLUGINDIR + '/madgraph/iolibs/template_files/'
    # add template file/ linking only needed in the madevent mode and not in standalone
    from_template = dict(PLUGIN_ProcessExporter.from_template)
    from_template['SubProcesses'] = from_template['SubProcesses'] + [s+'gpu/fbridge_common.inc',
                                      s+'gpu/counters.cc',
                                      s+'gpu/ompnumthreads.cc']

    to_link_in_P = PLUGIN_ProcessExporter.to_link_in_P + ['fbridge_common.inc', 'counters.cc','ompnumthreads.cc']

#------------------------------------------------------------------------------------

class SIMD_ProcessExporter(PLUGIN_ProcessExporter_MadEvent):

    # Default class for the run_card to use
    run_card_class = launch_plugin.CPPRunCard

    def change_output_args(args, cmd):
        """ """
        #cmd._export_format = "madevent_forplugin"
        cmd._export_format = 'madevent'
        cmd._export_plugin = FortranExporterBridge
        args.append('--hel_recycling=False')
        args.append('--me_exporter=standalone_simd')
        if 'vector_size' not in ''.join(args):
            args.append('--vector_size=16')
        if 'nb_wrap' not in ''.join(args):
            args.append('--nb_wrap=1')
        return args

class FortranExporterBridge(export_v4.ProcessExporterFortranMEGroup):
    _file_path = export_v4._file_path

    def write_auto_dsig_file(self, writer, matrix_element, proc_id = ""):
        replace_dict,context = super().write_auto_dsig_file(False, matrix_element, proc_id)
        replace_dict['additional_header'] = """
      INTEGER IEXT

      INTEGER                    ISUM_HEL
      LOGICAL                    MULTI_CHANNEL
      COMMON/TO_MATRIX/ISUM_HEL, MULTI_CHANNEL

      LOGICAL FIRST_CHID
      SAVE FIRST_CHID
      DATA FIRST_CHID/.TRUE./

#ifdef MG5AMC_MEEXPORTER_CUDACPP
      INCLUDE 'coupl.inc' ! for ALL_G
      INCLUDE 'fbridge.inc'
      INCLUDE 'fbridge_common.inc'
      INCLUDE 'genps.inc'
      INCLUDE 'run.inc'
      DOUBLE PRECISION OUT2(VECSIZE_MEMMAX)
      INTEGER SELECTED_HEL2(VECSIZE_MEMMAX)
      INTEGER SELECTED_COL2(VECSIZE_MEMMAX)
      DOUBLE PRECISION CBYF1
      INTEGER*4 NGOODHEL, NTOTHEL

      INTEGER*4 NWARNINGS
      SAVE NWARNINGS
      DATA NWARNINGS/0/

      LOGICAL FIRST
      SAVE FIRST
      DATA FIRST/.TRUE./"""
        replace_dict['OMP_LIB'] = ''
        replace_dict['OMP_PREFIX'] = """IF( FBRIDGE_MODE .LE. 0 ) THEN ! (FortranOnly=0 or BothQuiet=-1 or BothDebug=-2)
#endif
CALL COUNTERS_SMATRIX1MULTI_START( -1, VECSIZE_USED )  ! fortranMEs=-1"""
        replace_dict["OMP_POSTFIX"] = open(pjoin(PLUGINDIR,'madgraph','iolibs','template_files','gpu','smatrix_multi.f')).read().split('\n',4)[4] # AV skip 4 copyright lines
        if writer:
            file = open(pjoin(self._file_path, 'iolibs/template_files/auto_dsig_v4.inc')).read()
            file = file % replace_dict
            # Write the file
            writer.writelines(file, context=context)
        else:
            return replace_dict, context

    def write_driver(self, writer, *args, **kwargs):
        """Write the SubProcess/driver.f file with additions from CUDACPP"""
        replace_dict = super().write_driver(False, *args, **kwargs)

        # Additions from CUDACPP plugin (after patch)
        replace_dict['DRIVER_EXTRA_HEADER'] += """
      character*255 env_name, env_value
      integer env_length, env_status

#ifdef MG5AMC_MEEXPORTER_CUDACPP
      INCLUDE 'fbridge.inc'
c     INCLUDE 'fbridge_common.inc'
#endif
      INCLUDE 'fbridge_common.inc'
"""

        replace_dict['DRIVER_EXTRA_INITIALISE'] += """
#ifdef _OPENMP
      CALL OMPNUMTHREADS_NOT_SET_MEANS_ONE_THREAD()
#endif
      CALL COUNTERS_INITIALISE()

#ifdef MG5AMC_MEEXPORTER_CUDACPP
      fbridge_mode = 1 ! CppOnly=1, default for CUDACPP
#else
      fbridge_mode = 0 ! FortranOnly=0, default for FORTRAN
#endif
      env_name = 'CUDACPP_RUNTIME_FBRIDGEMODE'
      call get_environment_variable(env_name, env_value, env_length, env_status)
      if( env_status.eq.0 ) then
        write(*,*) 'Found environment variable "', trim(env_name), '" with value "', trim(env_value), '"'
        read(env_value,'(I255)') FBRIDGE_MODE ! see https://gcc.gnu.org/onlinedocs/gfortran/ICHAR.html
        write(*,*) 'FBRIDGE_MODE (from env) = ', FBRIDGE_MODE
      else if( env_status.eq.1 ) then ! 1 = not defined
        write(*,*) 'FBRIDGE_MODE (default) = ', FBRIDGE_MODE
      else ! -1 = too long for env_value, 2 = not supported by O/S
        write(*,*) 'ERROR! get_environment_variable failed for "', trim(env_name), '"'
        STOP
      endif
#ifndef MG5AMC_MEEXPORTER_CUDACPP
      if( fbridge_mode.ne.0 ) then
        write(*,*) 'ERROR! Invalid fbridge_mode (in FORTRAN backend mode) = ', fbridge_mode
        STOP
      endif
#endif

      env_name = 'CUDACPP_RUNTIME_VECSIZEUSED'
      call get_environment_variable(env_name, env_value, env_length, env_status)
      if( env_status.eq.0 ) then
        write(*,*) 'Found environment variable "', trim(env_name), '" with value "', trim(env_value), '"'
        read(env_value,'(I255)') VECSIZE_USED ! see https://gcc.gnu.org/onlinedocs/gfortran/ICHAR.html
        write(*,*) 'VECSIZE_USED (from env) = ', VECSIZE_USED
      else if( env_status.eq.1 ) then ! 1 = not defined
        write(*,*) 'VECSIZE_USED (default) = ', VECSIZE_USED
      else ! -1 = too long for env_value, 2 = not supported by O/S
        write(*,*) 'ERROR! get_environment_variable failed for "', trim(env_name), '"'
        STOP
      endif
      if( VECSIZE_USED.gt.VECSIZE_MEMMAX .or. VECSIZE_USED.le.0 ) then
        write(*,*) 'ERROR! Invalid VECSIZE_USED = ', VECSIZE_USED
        STOP
      endif

#ifdef MG5AMC_MEEXPORTER_CUDACPP
      CALL FBRIDGECREATE(FBRIDGE_PBRIDGE, VECSIZE_USED, NEXTERNAL, 4) ! this must be at the beginning as it initialises the CUDA device
      FBRIDGE_NCBYF1 = 0
      FBRIDGE_CBYF1SUM = 0
      FBRIDGE_CBYF1SUM2 = 0
      FBRIDGE_CBYF1MAX = -1D100
      FBRIDGE_CBYF1MIN = 1D100
#endif
"""

        replace_dict['DRIVER_EXTRA_FINALISE'] += """
#ifdef MG5AMC_MEEXPORTER_CUDACPP
      CALL FBRIDGEDELETE(FBRIDGE_PBRIDGE) ! this must be at the end as it shuts down the CUDA device
      IF( FBRIDGE_MODE .LE. -1 ) THEN ! (BothQuiet=-1 or BothDebug=-2)
        WRITE(*,'(a,f10.8,a,e8.2)')
     &    ' [MERATIOS] ME ratio CudaCpp/Fortran: MIN = ',
     &    FBRIDGE_CBYF1MIN + 1, ' = 1 - ', -FBRIDGE_CBYF1MIN
        WRITE(*,'(a,f10.8,a,e8.2)')
     &    ' [MERATIOS] ME ratio CudaCpp/Fortran: MAX = ',
     &    FBRIDGE_CBYF1MAX + 1, ' = 1 + ', FBRIDGE_CBYF1MAX
        WRITE(*,'(a,i6)')
     &    ' [MERATIOS] ME ratio CudaCpp/Fortran: NENTRIES = ',
     &    FBRIDGE_NCBYF1
c        WRITE(*,'(a,e8.2)')
c    &    ' [MERATIOS] ME ratio CudaCpp/Fortran - 1: AVG = ',
c    &    FBRIDGE_CBYF1SUM / FBRIDGE_NCBYF1
c       WRITE(*,'(a,e8.2)')
c    &    ' [MERATIOS] ME ratio CudaCpp/Fortran - 1: STD = ',
c    &    SQRT( FBRIDGE_CBYF1SUM2 / FBRIDGE_NCBYF1 ) ! ~standard deviation
        WRITE(*,'(a,e8.2,a,e8.2)')
     &    ' [MERATIOS] ME ratio CudaCpp/Fortran - 1: AVG = ',
     &    FBRIDGE_CBYF1SUM / FBRIDGE_NCBYF1, ' +- ',
     &    SQRT( FBRIDGE_CBYF1SUM2 ) / FBRIDGE_NCBYF1 ! ~standard error
      ENDIF
#endif
      CALL COUNTERS_FINALISE()
"""

        if writer:
            text = open(pjoin(self._file_path,'iolibs','template_files','madevent_driver.f')).read() % replace_dict
            writer.write(text)
            return True
        return replace_dict
#------------------------------------------------------------------------------------

class GPU_ProcessExporter(PLUGIN_ProcessExporter_MadEvent):

    # Default class for the run_card to use
    run_card_class = launch_plugin.GPURunCard

    def change_output_args(args, cmd):
        """ """
        cmd._export_format = 'madevent'
        cmd._export_plugin = FortranExporterBridge

        args.append('--hel_recycling=False')
        args.append('--me_exporter=standalone_cuda')
        if 'vector_size' not in ''.join(args):
            args.append('--vector_size=32')
        if 'nb_wrap' not in ''.join(args):
            args.append('--nb_wrap=512')
        return args

    def finalize(self, matrix_element, cmdhistory, MG5options, outputflag):
        misc.sprint("enter dedicated function")
        out = super().finalize(matrix_element, cmdhistory, MG5options, outputflag)
        # OM change RunCard class to have default for GPU
        text = open(pjoin(self.dir_path, 'bin', 'internal', 'launch_plugin.py'), 'r').read()
        text = text.replace('RunCard = CPPRunCard', 'RunCard = GPURunCard')
        open(pjoin(self.dir_path, 'bin', 'internal', 'launch_plugin.py'), 'w').write(text)
        return out

#------------------------------------------------------------------------------------
