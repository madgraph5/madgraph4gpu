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

    ###aloha_writer = 'cudac' #this was the default mode assigned to GPU 
    aloha_writer = ALOHAWriterForGPU # this is equivalent to the above line but allow to edit it obviously
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

    def write_aloha_routines(self):
        """Generate the hel_amps_model.h and hel_amps_model.cc files, which
        have the complete set of generalized Helas routines for the model"""
        
        if not os.path.isdir(os.path.join(self.dir_path, self.include_dir)):
            os.makedirs(os.path.join(self.dir_path, self.include_dir))
        if not os.path.isdir(os.path.join(self.dir_path, self.cc_file_dir)):
            os.makedirs(os.path.join(self.dir_path, self.cc_file_dir))

        model_h_file = os.path.join(self.dir_path, self.include_dir,
                                    'HelAmps_%s.h' % self.model_name)
        model_cc_file = os.path.join(self.dir_path, self.cc_file_dir,
                                     'HelAmps_%s.%s' % (self.model_name, self.cc_ext))

        replace_dict = {}

        replace_dict['output_name'] = self.output_name
        replace_dict['info_lines'] = export_cpp.get_mg5_info_lines()
        replace_dict['namespace'] = self.namespace
        replace_dict['model_name'] = self.model_name

        # Read in the template .h and .cc files, stripped of compiler
        # commands and namespaces
        template_h_files = self.read_aloha_template_files(ext = 'h')
        template_cc_files = self.read_aloha_template_files(ext = 'cc')

        import aloha.create_aloha as create_aloha
        aloha_model = create_aloha.AbstractALOHAModel(self.model.get('name'),
                                                      explicit_combine=True)
        aloha_model.add_Lorentz_object(self.model.get('lorentz'))
        
        if self.wanted_lorentz:
            aloha_model.compute_subset(self.wanted_lorentz)
        else:
            aloha_model.compute_all(save=False, custom_propa=True)
            
        for abstracthelas in dict(aloha_model).values():
            h_rout, cc_rout = abstracthelas.write(output_dir=None, 
                                                  language=self.aloha_writer, 
                                                  mode='no_include')

            template_h_files.append(h_rout)
            template_cc_files.append(cc_rout)
            
            #aloha_writer = aloha_writers.ALOHAWriterForCPP(abstracthelas,
            #                                               self.dir_path)
            #header = aloha_writer.define_header()
            #template_h_files.append(self.write_function_declaration(\
            #                             aloha_writer, header))
            #template_cc_files.append(self.write_function_definition(\
            #                              aloha_writer, header))

        replace_dict['function_declarations'] = '\n'.join(template_h_files)
        replace_dict['function_definitions'] = '\n'.join(template_cc_files)

        file_h = self.read_template_file(self.aloha_template_h) % replace_dict
        file_cc = self.read_template_file(self.aloha_template_cc) % replace_dict

        # Write the files
        import madgraph.iolibs.file_writers as writers
        ###writers.CPPWriter(model_h_file).writelines(file_h) # WITH FORMATTING
        ###writers.CPPWriter(model_cc_file).writelines(file_cc) # WITH FORMATTING
        writers.FileWriter(model_h_file).writelines(file_h) # WITHOUT FORMATTING
        writers.FileWriter(model_cc_file).writelines(file_cc) # WITHOUT FORMATTING

        import logging
        logger = logging.getLogger('madgraph.PLUGIN.CUDACPP_SA_OUTPUT.model_handling')
        logger.info("Created files %s and %s in directory" \
                    % (os.path.split(model_h_file)[-1],
                       os.path.split(model_cc_file)[-1]))
        logger.info("%s and %s" % \
                    (os.path.split(model_h_file)[0],
                     os.path.split(model_cc_file)[0]))

    def read_aloha_template_files(self, ext):
        """Read all ALOHA template files with extension ext, strip them of
        compiler options and namespace options, and return in a list"""
        # Use the plugin's path (for helas_h/cc)
        ###path = pjoin(MG5DIR, 'aloha','template_files')
        PLUGINDIR = os.path.dirname( __file__ )
        path = pjoin(PLUGINDIR, 'aloha', 'template_files')
        out = []
        if ext == 'h':
            out.append(open(pjoin(path, self.helas_h)).read())
        else:
            out.append(open(pjoin(path, self.helas_cc)).read())
        return out

    #===============================================================================
    # Global helper methods
    #===============================================================================
    @classmethod
    def read_template_file(cls, filename, classpath=False):
        """Open a template file and return the contents."""
        # Use the plugin's OneProcessExporterGPU template_path and __template_path (for aloha_template_h/cc)
        return OneProcessExporterGPU.read_template_file(filename, classpath)


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
    PLUGINDIR = os.path.dirname( __file__ )
    template_path = os.path.join( PLUGINDIR, 'madgraph', 'iolibs', 'template_files' )
    __template_path = os.path.join( PLUGINDIR, 'madgraph', 'iolibs', 'template_files' )
    process_template_h = 'gpu/process_h.inc'
    process_template_cc = 'gpu/process_cc.inc'
    process_class_template = 'gpu/process_class.inc'
    process_definition_template = 'gpu/process_function_definitions.inc'
    process_wavefunction_template = 'cpp_process_wavefunctions.inc'
    process_sigmaKin_function_template = 'gpu/process_sigmaKin_function.inc'
    single_process_template = 'gpu/process_matrix.inc'
    cc_ext = 'cu'
