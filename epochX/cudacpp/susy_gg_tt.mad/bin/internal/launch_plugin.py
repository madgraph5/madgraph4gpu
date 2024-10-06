# Copyright (C) 2020-2024 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: O. Mattelaer (Aug 2023) for the MG5aMC CUDACPP plugin.
# Further modified by: O. Mattelaer, A. Valassi (2024) for the MG5aMC CUDACPP plugin.

import logging
import os
import subprocess
pjoin = os.path.join
logger = logging.getLogger('cmdprint') # for stdout

try:
    import madgraph
except ImportError:
    import internal.madevent_interface as madevent_interface
    import internal.misc as misc
    import internal.extended_cmd as extended_cmd
    import internal.banner as banner_mod
    import internal.common_run_interface as common_run_interface
else:
    import madgraph.interface.madevent_interface as madevent_interface
    import madgraph.various.misc as misc
    import madgraph.interface.extended_cmd as extended_cmd
    import madgraph.various.banner as banner_mod
    import madgraph.interface.common_run_interface as common_run_interface

class CPPMEInterface(madevent_interface.MadEventCmdShell):
    def compile(self, *args, **opts):
        """ """
        import multiprocessing
        if not self.options['nb_core'] or self.options['nb_core'] == 'None':
            self.options['nb_core'] = multiprocessing.cpu_count()    
        if 'cwd' in opts and os.path.basename(opts['cwd']) == 'Source':
            path = pjoin(opts['cwd'], 'make_opts')
            common_run_interface.CommonRunCmd.update_make_opts_full(path,
                {'FPTYPE': self.run_card['cudacpp_fptype'],
                 'HELINL': self.run_card['cudacpp_helinl'],
                 'HRDCOD': self.run_card['cudacpp_hrdcod'] })
            misc.sprint('FPTYPE, HELINL, HRDCOD checked')
        cudacpp_supported_backends = [ 'fortran', 'cuda', 'hip', 'cpp', 'cppnone', 'cppsse4', 'cppavx2', 'cpp512y', 'cpp512z', 'cppauto' ]
        if args and args[0][0] == 'madevent' and hasattr(self, 'run_card'):            
            if self.run_card['cudacpp_bldall'] == True: # pre-build all backends #945
                logger.info("Pre-building madevent in madevent_interface.py with ALL matrix elements")
                args[0][0] = 'bldall'
                misc.compile(nb_core=self.options['nb_core'], *args, **opts)
            cudacpp_backend = self.run_card['cudacpp_backend'].lower() # the default value is defined in launch_plugin.py
            logger.info("Building madevent in madevent_interface.py with '%s' matrix elements"%cudacpp_backend)
            if cudacpp_backend in cudacpp_supported_backends :
                args[0][0] = 'madevent_' + cudacpp_backend + '_link'
            else:
                raise Exception( "Invalid cudacpp_backend='%s': supported backends are %s"%supported_backends )
            return misc.compile(nb_core=self.options['nb_core'], *args, **opts)
        else:
            return misc.compile(nb_core=self.options['nb_core'], *args, **opts)

# CUDACPP runcard block ------------------------------------------------------------------------------------
template_on = \
"""#***********************************************************************
# SIMD/GPU configuration for the CUDACPP plugin
#************************************************************************
 %(cudacpp_backend)s = cudacpp_backend ! CUDACPP backend: fortran, cuda, hip, cpp (DEFAULT), cppnone, cppsse4, cppavx2, cpp512y, cpp512z, cppauto
#*** WARNING! The following cudacpp runcards are experimental! Users should normally change only the cudacpp_backend card ***
 %(cudacpp_fptype)s = cudacpp_fptype ! CUDACPP floating point precision: f (single), d (double), m (mixed, DEFAULT: double for amplitudes, single for colors)
 %(cudacpp_hrdcod)s = cudacpp_hrdcod ! CUDACPP parameter hardcoding: 0 (DEFAULT, parameters not hardcoded: read param_card.dat at runtime), 1 (hardcoded parameters)
 %(cudacpp_helinl)s = cudacpp_helinl ! CUDACPP helicity amplitude inlining: 0 (DEFAULT, ordinary inlining of templates), 1 (aggressive inlining with 'always inline')
 %(cudacpp_bldall)s = cudacpp_bldall ! CUDACPP build all available backends in separate build directories: False, True
"""
template_off = ''
plugin_block = banner_mod.RunBlock('simd', template_on=template_on, template_off=template_off)

class CPPRunCard(banner_mod.RunCardLO):
    blocks = banner_mod.RunCardLO.blocks + [plugin_block]

    def reset_simd(self, old_value, new_value, name):
        if not hasattr(self, 'path'):
            raise Exception('INTERNAL ERROR! CPPRunCard instance has no attribute path') # now ok after fixing #790
        if name == "vector_size" and new_value <= int(old_value):
            # code can handle the new size -> do not recompile
            return

        # ok need to force recompilation of the cpp part
        Sourcedir = pjoin(os.path.dirname(os.path.dirname(self.path)), 'Source')
        subprocess.call(['make', 'cleanall'], cwd=Sourcedir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    def reset_makeopts(self, old_value, new_value, name):
        if not hasattr(self, 'path'):
            raise Exception
        if name == 'cudacpp_fptype':
            common_run_interface.CommonRunCmd.update_make_opts_full({'FPTYPE': new_value})
        elif name == 'cudacpp_hrdcod':
            raise Exception('Cannot change cudacpp_hrdcod')
        elif name == 'cudacpp_helinl':
            raise Exception('Cannot change cudacpp_helinl')
        elif name == 'cudacpp_bldall':
            raise Exception('Cannot change cudacpp_bldall')
        else:
            raise Exception
        Sourcedir = pjoin(os.path.dirname(os.path.dirname(self.path)), 'Source')
        subprocess.call(['make', 'cleanall'], cwd=Sourcedir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    def plugin_input(self, finput):
        return

    def default_setup(self):
        super().default_setup()
        cudacpp_supported_backends = [ 'fortran', 'cuda', 'hip', 'cpp', 'cppnone', 'cppsse4', 'cppavx2', 'cpp512y', 'cpp512z', 'cppauto' ]
        self.add_param('cudacpp_backend', 'cpp',
                       include=False, # AV: 'include=True' would add "CUDACPP_BACKEND = 'cpp'" to run_card.inc
                       hidden=False, # AV: keep cudacpp_backend in runcard template and keep 'hidden='False'
                       allowed=cudacpp_supported_backends)
        self.add_param('cudacpp_fptype', 'm',
                       include=False, # AV: 'include=True' would add "CUDACPP_FPTYPE = 'm'" to run_card.inc (if fct_mod is removed, else codegen fails)
                       hidden=False, # AV: add cudacpp_backend to runcard template and keep 'hidden='False'
                       fct_mod=(self.reset_makeopts,(),{}), # AV: I assume this forces a 'make cleanavx' if FPTYPE changes?
                       allowed=['m','d','f']
                       )
        self.add_param('cudacpp_helinl', '0',
                       include=False,  # AV: no need to add this parameter to run_card.inc
                       hidden=False, # AV: add cudacpp_helinl to runcard template and keep 'hidden='False'
                       fct_mod=(self.reset_makeopts,(),{}), # AV: I assume this raises an exception if cudacpp_helinl changes?
                       allowed=['0','1']
                       )
        self.add_param('cudacpp_hrdcod', '0',
                       include=False,  # AV: no need to add this parameter to run_card.inc
                       hidden=False, # AV: add cudacpp_hrdcod to runcard template and keep 'hidden='False'
                       fct_mod=(self.reset_makeopts,(),{}), # AV: I assume this raises an exception if cudacpp_hrdcod changes?
                       allowed=['0','1']
                       )
        self.add_param('cudacpp_bldall', False,
                       include=False, # AV: no need to add this parameter to run_card.inc
                       hidden=False, # AV: add cudacpp_bldall to runcard template and keep 'hidden='False'
                       fct_mod=(self.reset_makeopts,(),{}), # AV: I assume this raises an exception if cudacpp_bldall changes?
                       )
        self['vector_size'] = 16 # already setup in default class (just change value)
        self['aloha_flag'] = '--fast-math'
        self['matrix_flag'] = '-O3'
        self['limhel'] = 0
        self.display_block.append('simd')
        self.display_block.append('psoptim')

    # OM/AV - overload the default version in banner.py
    def write_one_include_file(self, output_dir, incname, output_file=None):
        """write one include file at the time"""
        if incname == "vector.inc":
            if 'vector_size' not in self.user_set and 'wrap_size' not in self.user_set: return
            if output_file is None: vectorinc=pjoin(output_dir,incname)
            else: vectorinc=output_file
            with open(vectorinc+'.new','w') as fileout:
                with open(vectorinc) as filein:
                    for line in filein:
                        if line.startswith('C'): fileout.write(line)
            super().write_one_include_file(output_dir, incname, output_file)
            with open(vectorinc+'.new','a') as fileout:
                with open(vectorinc) as filein:
                    for line in filein:
                        if not line.startswith('\n'): fileout.write(line)
            os.replace(vectorinc+'.new',vectorinc)
        else:
            super().write_one_include_file(output_dir, incname, output_file)

    def check_validity(self):
        """ensure that PLUGIN information are consistent"""
        super().check_validity()
        if self['SDE_strategy'] != 1:
            logger.warning('SDE_strategy different of 1 is not supported with SMD/GPU mode')
            self['sde_strategy'] = 1
        if self['hel_recycling']:
            self['hel_recycling'] = False

class GPURunCard(CPPRunCard):
    def default_setup(self):
        super().default_setup()
        # change default value:
        self['cudacpp_backend'] = 'cuda'
        self['vector_size'] = 16384 # already setup in default class (just change value)

MEINTERFACE = CPPMEInterface
RunCard = CPPRunCard
