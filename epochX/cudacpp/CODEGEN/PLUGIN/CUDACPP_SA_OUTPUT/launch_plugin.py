# Copyright (C) 2020-2025 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: O. Mattelaer (Aug 2023) for the MG5aMC CUDACPP plugin.
# Further modified by: O. Mattelaer, A. Valassi, Z. Wettersten (2024-2025) for the MG5aMC CUDACPP plugin.

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
    import internal.files as files
else:
    import madgraph.interface.madevent_interface as madevent_interface
    import madgraph.various.misc as misc
    import madgraph.interface.extended_cmd as extended_cmd
    import madgraph.various.banner as banner_mod
    import madgraph.interface.common_run_interface as common_run_interface
    import madgraph.iolibs.files as files

class CPPMEInterface(madevent_interface.MadEventCmdShell):
    def compile(self, *args, **opts):
        """ """
        import multiprocessing
        if not self.options['nb_core'] or self.options['nb_core'] == 'None':
            self.options['nb_core'] = multiprocessing.cpu_count()    
        if 'cwd' in opts and os.path.basename(opts['cwd']) == 'Source':
            path = pjoin(opts['cwd'], 'make_opts')
            common_run_interface.CommonRunCmd.update_make_opts_full(path,
                {'override FPTYPE': self.run_card['floating_type'] })
            misc.sprint('FPTYPE checked')
        cudacpp_supported_backends = [ 'fortran', 'cuda', 'hip', 'cpp', 'cppnone', 'cppsse4', 'cppavx2', 'cpp512y', 'cpp512z', 'cppauto' ]
        if args and args[0][0] == 'madevent' and hasattr(self, 'run_card'):            
            cudacpp_backend = self.run_card['cudacpp_backend'].lower() # the default value is defined in launch_plugin.py
            if cudacpp_backend in ['cpp', 'cppauto']:
                backend_log = pjoin(opts["cwd"], ".resolved-backend")
                # try to remove old file if present
                try:
                    os.remove(backend_log)
                except FileNotFoundError:
                    pass
                misc.compile(["-f", "cudacpp.mk", f"BACKEND=cppauto", f"BACKEND_LOG={backend_log}", "detect-backend"], **opts)
                try:
                    with open(backend_log, "r") as f:
                        resolved_backend = f.read().strip()
                    logger.info(f"Backend '{cudacpp_backend}' resolved as '{resolved_backend}'")
                    cudacpp_backend = resolved_backend
                except FileNotFoundError:
                    raise RuntimeError("Could not resolve cudacpp_backend=cppauto|cpp; ensure Makefile detection runs properly.")
            logger.info(f"Building madevent in madevent_interface.py with '{cudacpp_backend}' matrix elements")
            if cudacpp_backend in cudacpp_supported_backends :
                args[0][0] = 'madevent_' + cudacpp_backend + '_link'
            else:
                raise Exception(f"Invalid cudacpp_backend='{cudacpp_backend}': supported backends are [ '" + "', '".join(cudacpp_supported_backends) + "' ]")
            return misc.compile(nb_core=self.options['nb_core'], *args, **opts)
        else:
            return misc.compile(nb_core=self.options['nb_core'], *args, **opts)

    def do_generate_events(self, *args, **kwargs):
        cudacpp_version = os.path.join(self.me_dir, "CUDACPP_VERSION.txt")
        if os.path.exists(cudacpp_version):
            with open(cudacpp_version, "r") as f:
                lines = f.readlines()
            logger.info("=================================================")
            for line in lines:
                logger.info(line.strip())
            logger.info("=================================================")
        return super().do_generate_events(*args, **kwargs)

    def do_create_gridpack(self, *args, **kwargs):
        """Overload to embed the CUDACPP_VERSION.txt banner into the gridpack.
        The banner is not printed here (this runs at gridpack *creation* time); instead
        it is packaged inside the tarball and printed at *run* time by run.sh, since a
        gridpack run uses GridPackCmd and never goes through this plugin interface."""
        self.embed_cudacpp_version_in_gridpack()
        return super().do_create_gridpack(*args, **kwargs)

    # DM - make the CUDACPP_VERSION.txt banner available (and printed) inside gridpacks
    def embed_cudacpp_version_in_gridpack(self):
        """Prepare the process directory so that make_gridpack packages the CUDACPP banner
        and run.sh prints it at runtime. Only files copied into the tarball are touched:
          - copy CUDACPP_VERSION.txt into bin/internal/ (bin/ is packaged into madevent/)
          - patch bin/internal/Gridpack/run.sh (the template packaged as ./run.sh) to cat it
        This must run *before* super().do_create_gridpack(), which invokes make_gridpack."""
        version_src = pjoin(self.me_dir, 'CUDACPP_VERSION.txt')
        if not os.path.exists(version_src):
            logger.warning('CUDACPP_VERSION.txt not found in %s: the gridpack will not print '
                           'the CUDACPP version banner' % self.me_dir)
            return
        # 1) copy the banner into a directory that make_gridpack moves into madevent/
        #    (make_gridpack packages 'bin' but not the process-root CUDACPP_VERSION.txt)
        files.cp(version_src, pjoin(self.me_dir, 'bin', 'internal', 'CUDACPP_VERSION.txt'))
        # 2) patch the gridpack run.sh so it prints the banner at runtime (idempotent)
        runsh = pjoin(self.me_dir, 'bin', 'internal', 'Gridpack', 'run.sh')
        if not os.path.exists(runsh):
            logger.warning('%s not found: the gridpack will not print the CUDACPP version banner'
                           % runsh)
            return
        marker = '# CUDACPP version banner'
        with open(runsh) as fsock:
            content = fsock.read()
        if marker in content:
            return # already patched (e.g. create_gridpack called more than once)
        # DIR (=./madevent when the gridpack is unpacked) is defined just above this anchor
        anchor = '# For Linux'
        banner_block = (
            '%s (printed by the MG5aMC CUDACPP plugin)\n'
            'if [ -f "${DIR}/bin/internal/CUDACPP_VERSION.txt" ]; then\n'
            '    echo "================================================="\n'
            '    cat "${DIR}/bin/internal/CUDACPP_VERSION.txt"\n'
            '    echo "================================================="\n'
            'fi\n\n' % marker
        )
        if anchor not in content:
            logger.warning('Could not find the expected anchor in %s: the gridpack will not '
                           'print the CUDACPP version banner' % runsh)
            return
        content = content.replace(anchor, banner_block + anchor, 1)
        with open(runsh, 'w') as fsock:
            fsock.write(content)
        logger.info('Patched %s to print the CUDACPP version banner at gridpack runtime' % runsh)

# Phase-Space Optimization ------------------------------------------------------------------------------------
template_on = \
"""#***********************************************************************
# SIMD/GPU configuration for the CUDACPP plugin
#************************************************************************
 %(cudacpp_backend)s = cudacpp_backend ! CUDACPP backend: fortran, cuda, hip, cpp, cppnone, cppsse4, cppavx2, cpp512y, cpp512z, cppauto
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
        if name == 'floating_type':
            common_run_interface.CommonRunCmd.update_make_opts_full({'override FPTYPE': new_value})
        else:
            raise Exception
        Sourcedir = pjoin(os.path.dirname(os.path.dirname(self.path)), 'Source')
        subprocess.call(['make', 'cleanall'], cwd=Sourcedir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    def default_setup(self):
        super().default_setup()
        self.add_param('floating_type', 'm', include=False, hidden=True,
                       fct_mod=(self.reset_makeopts,(),{}),
                       allowed=['m','d','f'],
                       comment='floating point precision: f (single), d (double), m (mixed: double for amplitudes, single for colors)'
                       )
        cudacpp_supported_backends = [ 'fortran', 'cuda', 'hip', 'cpp', 'cppnone', 'cppsse4', 'cppavx2', 'cpp512y', 'cpp512z', 'cppauto' ]
        self.add_param('cudacpp_backend', 'cpp', include=False, hidden=False,
                       allowed=cudacpp_supported_backends)
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
        self['vector_size'] = 32 # ZW: default to 32, might want to change to 64 to utilise AMD GPUs better as well # 16384 # already setup in default class (just change value)
        self['nb_warp'] = 512 # number of warps per kernel call, for now setting to 16 384 / vector_size

MEINTERFACE = CPPMEInterface
RunCard = CPPRunCard
