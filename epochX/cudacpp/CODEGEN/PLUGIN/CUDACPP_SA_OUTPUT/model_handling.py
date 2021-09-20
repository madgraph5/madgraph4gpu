import madgraph.iolibs.export_cpp as export_cpp
import aloha.aloha_writers as aloha_writers

import os
pjoin = os.path.join

class ALOHAWriterForGPU(aloha_writers.ALOHAWriterForGPU):
    
    extension = '.cu'
    prefix ='__device__'
    realoperator = '.real()'
    imagoperator = '.imag()'
    ci_definition = 'cxtype cI = cxtype(0., 1.);\n'
    
    type2def = {}    
    type2def['int'] = 'int '
    type2def['double'] = 'fptype '
    type2def['complex'] = 'cxtype '
    type2def['pointer_vertex'] = '*' # using complex<double> * vertex)
    type2def['pointer_coup'] = ''


class  UFOModelConverterGPU(export_cpp.UFOModelConverterGPU):

    #aloha_writer = 'cudac' #this was the default mode assigned to GPU 
    aloha_writer = ALOHAWriterForGPU # this is equivalent to the above line but allow to edit it obviously.
    cc_ext = 'cu'
    # Template files to use
    #include_dir = '.'
    #c_file_dir = '.'
    #param_template_h = 'cpp_model_parameters_h.inc'
    #param_template_cc = 'cpp_model_parameters_cc.inc'
    aloha_template_h = pjoin('gpu','cpp_hel_amps_h.inc')
    aloha_template_cc = pjoin('gpu','cpp_hel_amps_cc.inc')
    helas_h = pjoin('gpu', 'helas.h')
    helas_cc = pjoin('gpu', 'helas.cu')


import madgraph.iolibs.helas_call_writers as helas_call_writers
    
class GPUFOHelasCallWriter(helas_call_writers.GPUFOHelasCallWriter):

    def format_coupling(self, call):
        """Format the coupling so any minus signs are put in front"""
        return super().format_coupling(call)
        

    def get_external(self,wf, argument):
        """ formatting for ixxxx/ oxxxx /.... type of function (external ones) """
        return super().get_external(wf, argument)

    def generate_helas_call(self, argument):
        """Routine for automatic generation of C++ Helas calls
        according to just the spin structure of the interaction.

        First the call string is generated, using a dictionary to go
        from the spin state of the calling wavefunction and its
        mothers, or the mothers of the amplitude, to difenrentiate wich call is
        done.

        Then the call function is generated, as a lambda which fills
        the call string with the information of the calling
        wavefunction or amplitude. The call has different structure,
        depending on the spin of the wavefunction and the number of
        mothers (multiplicity of the vertex). The mother
        wavefunctions, when entering the call, must be sorted in the
        correct way - this is done by the sorted_mothers routine.

        Finally the call function is stored in the relevant
        dictionary, in order to be able to reuse the function the next
        time a wavefunction with the same Lorentz structure is needed.
        """
        return super().generate_helas_call(argument)


class OneProcessExporterGPU(export_cpp.OneProcessExporterGPU):

        # Static variables (for inheritance)
    process_dir = '.'
    include_dir = '.'
    #template_path = os.path.join(_file_path, 'iolibs', 'template_files')
    #__template_path = os.path.join(_file_path, 'iolibs', 'template_files') 
    process_template_h = 'gpu/process_h.inc'
    process_template_cc = 'gpu/process_cc.inc'
    process_class_template = 'gpu/process_class.inc'
    process_definition_template = 'gpu/process_function_definitions.inc'
    process_wavefunction_template = 'cpp_process_wavefunctions.inc'
    process_sigmaKin_function_template = 'gpu/process_sigmaKin_function.inc'
    single_process_template = 'gpu/process_matrix.inc'
    cc_ext = 'cu'
