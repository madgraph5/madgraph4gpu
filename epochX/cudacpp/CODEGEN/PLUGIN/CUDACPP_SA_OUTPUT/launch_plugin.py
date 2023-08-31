
import madgraph.various.misc as misc
import madgraph.interface.extended_cmd as extended_cmd
import logging

logger = logging.getLogger('cmdprint') # for stdout

try:
    import madgraph
except:
    import internal.madevent_interface as madevent_interface
else:
    import madgraph.interface.madevent_interface as madevent_interface

class CPPMEInterface(madevent_interface.MadEventCmdShell):
    
    def compile(self, *args, **opts):
        """ """
        import multiprocessing
        if not self.options['nb_core'] or self.options['nb_core'] == 'None':
            self.options['nb_core'] = multiprocessing.cpu_count()
    
        if args and args[0][0] == 'madevent' and hasattr(self, 'run_card'):
            cudacpp_backend = self.run_card['cudacpp_backend'] # the default value is defined in banner.py
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
    
MEINTERFACE = CPPMEInterface
