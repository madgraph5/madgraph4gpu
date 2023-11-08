# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: O. Mattelaer (Aug 2023) for the MG5aMC CUDACPP plugin.
# Further modified by: O. Mattelaer, A. Valassi (2023) for the MG5aMC CUDACPP plugin.

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
else:
    import madgraph.interface.madevent_interface as madevent_interface
    import madgraph.various.misc as misc
    import madgraph.interface.extended_cmd as extended_cmd
    import madgraph.various.banner as banner_mod

class CPPMEInterface(madevent_interface.MadEventCmdShell):
    def compile(self, *args, **opts):
        """ """
        import multiprocessing
        if not self.options['nb_core'] or self.options['nb_core'] == 'None':
            self.options['nb_core'] = multiprocessing.cpu_count()
        if args and args[0][0] == 'madevent' and hasattr(self, 'run_card'):
            import pathlib
            import os
            pjoin = os.path.join
            cudacpp_backend = self.run_card['cudacpp_backend'].upper() # the default value is defined in banner.py
            logger.info("Building madevent in madevent_interface.py with '%s' matrix elements"%cudacpp_backend)
            if cudacpp_backend == 'FORTRAN':
                args[0][0] = 'madevent_fortran_link'
            elif cudacpp_backend == 'CPP':
                args[0][0] = 'madevent_cpp_link'
            elif cudacpp_backend == 'CUDA':
                args[0][0] = 'madevent_cuda_link'
            else:
                raise Exception("Invalid cudacpp_backend='%s': only 'FORTRAN', 'CPP', 'CUDA' are supported")
            return misc.compile(nb_core=self.options['nb_core'], *args, **opts)
        else:
            return misc.compile(nb_core=self.options['nb_core'], *args, **opts)

class CPPRunCard(banner_mod.RunCardLO):
    def reset_simd(self, old_value, new_value, name):
        if not hasattr(self, 'path'):
            raise Exception

        if name != 'vecsize_memmax':
            # this will be control by that value only
            return
        
        if name == "vecsize_memax" and new_value <= int(old_value):
            # code can handle the new size -> do not recompile
            return

        # ok need to force recompilation of the cpp part
        Sourcedir = pjoin(os.path.dirname(os.path.dirname(self.path)), 'Source')
        subprocess.call(['make', 'cleanavx'], cwd=Sourcedir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
    def plugin_input(self, finput):
        return

    def default_setup(self):
        super().default_setup()
        self.add_param('cudacpp_backend', 'CPP', include=False, hidden=False)

    def write_one_include_file(self, output_dir, incname, output_file=None):
        """write one include file at the time"""

        if incname == "vector.inc" and 'vector_size' not in self.user_set and\
           'wrap_size' not in self.user_set:
            return
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
        super(CPPRunCard, self).default_setup()
        self.add_param('cudacpp_backend', 'CUDA', include=False, hidden=False)

#class CUDACPPRunCard(CPPRunCard):
#    def default_setup(self):
#        super(CPPRunCard, self).default_setup()
#        self.add_param('cudacpp_backend', 'CPP', include=False, hidden=False)

MEINTERFACE = CPPMEInterface
RunCard = CPPRunCard
