import madgraph.iolibs.export_cpp as export_cpp
import aloha.aloha_writers as aloha_writers

import os
pjoin = os.path.join

import aloha
from six import StringIO
from collections import defaultdict

class ALOHAWriterForGPU(aloha_writers.ALOHAWriterForGPU):
    
    extension = '.cu'
    prefix ='__device__'
    realoperator = '.real()'
    imagoperator = '.imag()'
    ci_definition = 'cxtype cI = cxtype(0., 1.);\n'
    
    type2def = {}    
    type2def['int'] = 'int'
    type2def['double'] = 'fptype'
    type2def['complex'] = 'cxtype'
    type2def['pointer_vertex'] = '*' # using complex<double> * vertex)
    type2def['pointer_coup'] = ''

    def change_number_format(self, number):
        """Formating the number"""
        def isinteger(x):
            try:
                return int(x) == x
            except TypeError:
                return False
        if isinteger(number):
            if number == 1:
                out = 'one'
            elif number == -1:
                out = '-one'
            else:
                out = '%s.' % (str(int(number))) # This prints -1 as '-1.'
        elif isinstance(number, complex):
            if number.imag:
                if number.real:
                    out = '(%s + %s*cI)' % (self.change_number_format(number.real), \
                                    self.change_number_format(number.imag))
                else:
                    if number.imag == 1:
                        out = 'cI'
                    elif number.imag == -1:
                        out = '-cI'
                    else: 
                        out = '%s * cI' % self.change_number_format(number.imag)
            else:
                out = '%s' % (self.change_number_format(number.real))
        else:
            tmp = Fraction(str(number))
            tmp = tmp.limit_denominator(100)
            if not abs(tmp - number) / abs(tmp + number) < 1e-8:
                out = '%.9f' % (number)
            else:
                out = '%s./%s.' % (tmp.numerator, tmp.denominator)
        return out

    def get_header_txt(self, name=None, couplings=None,mode=''):
        """Define the Header of the fortran file. This include
            - function tag
            - definition of variable
        """
        if name is None:
            name = self.name
        if mode=='':
            mode = self.mode
        out = StringIO()
        # define the type of function and argument
        if not 'no_include' in mode:
            out.write('#include \"%s.h\"\n\n' % self.name)
        args = []
        for format, argname in self.define_argument_list(couplings):
            if format.startswith('list'):
                type = self.type2def[format[5:]]
                list_arg = '[]'
            else:
                type = self.type2def[format]
                list_arg = ''
            if argname.startswith('COUP'):
                point = self.type2def['pointer_coup']
                args.append('%s %s%s%s'% (type, point, argname, list_arg))
            else:
                args.append('%s %s%s'% (type, argname, list_arg))
        if not self.offshell:
            output = '%(doublec)s %(pointer_vertex)s vertex' % {
                'doublec':self.type2def['complex'],
                'pointer_vertex': self.type2def['pointer_vertex']}
            #self.declaration.add(('complex','vertex'))
        else:
            output = '%(doublec)s %(spin)s%(id)d[]' % {
                     'doublec': self.type2def['complex'],
                     'spin': self.particles[self.outgoing -1],
                     'id': self.outgoing}
            self.declaration.add(('list_complex', output))
        out.write('%(prefix)s void %(name)s(const %(args)s, %(output)s)' % \
                  {'prefix': self.prefix,
                      'output':output, 'name': name, 'args': ', const '.join(args)})
        if 'is_h' in mode:
            out.write(';\n')
        else:
            out.write('\n{\n')
        return out.getvalue() 

    def get_declaration_txt(self, add_i=True):
        """ Prototype for how to write the declaration of variable
            Include the symmetry line (entry FFV_2)
        """        
        out = StringIO()
        argument_var = [name for type,name in self.call_arg]
        # define the complex number CI = 0+1j
        if add_i:
            out.write('  ' + self.ci_definition)
        for type, name in self.declaration.tolist():
            if type.startswith('list'):
                type = type[5:]
                if name.startswith('P'):
                    size = 4
                elif not 'tmp' in name:
                    continue
                    #should be define in the header
                elif name[0] in ['F','V']:
                    if aloha.loop_mode:
                        size = 8
                    else:
                        size = 6
                elif name[0] == 'S':
                    if aloha.loop_mode:
                        size = 5
                    else:
                        size = 3
                elif name[0] in ['R','T']: 
                    if aloha.loop_mode:
                        size = 20
                    else:
                        size = 18
                out.write('  %s %s[%s];\n' % (self.type2def[type], name, size))
            elif (type, name) not in self.call_arg:
                out.write('  %s %s;\n' % (self.type2def[type], name))               
        return out.getvalue()

    def get_momenta_txt(self):
        """Define the Header of the fortran file. This include
            - momentum conservation
            - definition of the impulsion"""
        out = StringIO()
        # Define all the required momenta
        p = [] # a list for keeping track how to write the momentum
        signs = self.get_momentum_conservation_sign()
        for i,type in enumerate(self.particles):
            if self.declaration.is_used('OM%s' % (i+1)):
                out.write("    OM{0} = {1};\n    if (M{0} != {1})\n OM{0}={2}/(M{0}*M{0});\n".format( 
                         i+1, self.change_number_format(0), self.change_number_format(1)))
            if i+1 == self.outgoing:
                out_type = type
                out_size = self.type_to_size[type] 
                continue
            elif self.offshell:
                if len(p) != 0 : p.append(' ')
                p.append('{0} {1}{2}[%(i)s]'.format(signs[i],type,i+1,type))    
            if self.declaration.is_used('P%s' % (i+1)):
                self.get_one_momenta_def(i+1, out)
        # define the resulting momenta
        if self.offshell:
            energy_pos = out_size -2
            type = self.particles[self.outgoing-1]
            if aloha.loop_mode:
                size_p = 4
            else:
                size_p = 2
            for i in range(size_p):
                dict_energy = {'i':i}
                out.write('  %s%s[%s] = %s;\n' % (type,self.outgoing, i, 
                                             ''.join(p) % dict_energy))
            if self.declaration.is_used('P%s' % self.outgoing):
                self.get_one_momenta_def(self.outgoing, out)
        # Returning result
        return out.getvalue()

    def get_one_momenta_def(self, i, strfile):
        type = self.particles[i-1]
        if aloha.loop_mode:
            template ='  P%(i)d[%(j)d] = %(sign)s %(type)s%(i)d[%(nb)d];\n'
        else:
            template ='  P%(i)d[%(j)d] = %(sign)s %(type)s%(i)d[%(nb2)d]%(operator)s;\n'
        nb2 = 0
        for j in range(4):
            if not aloha.loop_mode:
                nb = j 
                if j == 0: 
                    assert not aloha.mp_precision 
                    operator = self.realoperator # not suppose to pass here in mp
                elif j == 1: 
                    nb2 += 1
                elif j == 2:
                    assert not aloha.mp_precision 
                    operator = self.imagoperator # not suppose to pass here in mp
                elif j ==3:
                    nb2 -= 1
            else:
                operator =''
                nb = j
                nb2 = j
            sign = self.get_P_sign(i) if self.get_P_sign(i) else '+'
            strfile.write(template % {'j':j,'type': type, 'i': i, 
                        'nb': nb, 'nb2': nb2, 'operator':operator,
                        'sign': sign})

    def define_expression(self):
        """Write the helicity amplitude in C++ format"""
        out = StringIO()
        if self.routine.contracted:
            keys = sorted(self.routine.contracted.keys())
            for name in keys:
                obj = self.routine.contracted[name]
                out.write('  %s = %s;\n' % (name, self.write_obj(obj)))
                self.declaration.add(('complex', name))
        for name, (fct, objs) in self.routine.fct.items():
            format = ' %s = %s;\n' % (name, self.get_fct_format(fct))
            out.write(format % ','.join([self.write_obj(obj) for obj in objs]))
        numerator = self.routine.expr
        if not 'Coup(1)' in self.routine.infostr:
            coup_name = 'COUP'
        else:
            coup_name = '%s' % self.change_number_format(1)
        if not self.offshell:
            if coup_name == 'COUP':
                mydict = {'num': self.write_obj(numerator.get_rep([0]))}
                for c in ['coup', 'vertex']:
                    if self.type2def['pointer_%s' %c] in ['*']:
                        mydict['pre_%s' %c] = '(*'
                        mydict['post_%s' %c] = ')'
                    else:
                        mydict['pre_%s' %c] = ''
                        mydict['post_%s'%c] = ''
                out.write('  %(pre_vertex)svertex%(post_vertex)s = %(pre_coup)sCOUP%(post_coup)s * %(num)s;\n' %\
                            mydict)
            else:
                mydict= {}
                if self.type2def['pointer_vertex'] in ['*']:
                    mydict['pre_vertex'] = '(*'
                    mydict['post_vertex'] = ')'
                else:
                    mydict['pre_vertex'] = ''
                    mydict['post_vertex'] = ''                 
                mydict['data'] = self.write_obj(numerator.get_rep([0]))
                out.write(' %(pre_vertex)svertex%(post_vertex)s = %(data)s;\n' % 
                          mydict)
        else:
            OffShellParticle = '%s%d' % (self.particles[self.offshell-1],\
                                                                  self.offshell)
            if 'L' not in self.tag:
                coeff = 'denom'
                mydict = {}
                if self.type2def['pointer_coup'] in ['*']:
                    mydict['pre_coup'] = '(*'
                    mydict['post_coup'] = ')'
                else:
                    mydict['pre_coup'] = ''
                    mydict['post_coup'] = ''
                mydict['coup'] = coup_name
                mydict['i'] = self.outgoing
                if not aloha.complex_mass:
                    if self.routine.denominator:
                        out.write('  denom = %(pre_coup)s%(coup)s%(post_coup)s / (%(denom)s)\n' % \
                                  mydict) 
                    else:
                        out.write('  denom = %(pre_coup)s%(coup)s%(post_coup)s / ((P%(i)s[0]*P%(i)s[0]) - (P%(i)s[1]*P%(i)s[1]) - (P%(i)s[2]*P%(i)s[2]) - (P%(i)s[3]*P%(i)s[3]) - M%(i)s*(M%(i)s-cI*W%(i)s));\n' % \
                                  mydict)
                else:
                    if self.routine.denominator:
                        raise Exception('modify denominator are not compatible with complex mass scheme')                
                    out.write('  denom = %(pre_coup)s%(coup)s%(post_coup)s / ((P%(i)s[0]*P%(i)s[0])-(P%(i)s[1]*P%(i)s[1])-(P%(i)s[2]*P%(i)s[2])-(P%(i)s[3]*P%(i)s[3]) - (M%(i)s*M%(i)s));\n' % \
                              mydict)
                self.declaration.add(('complex','denom'))
                if aloha.loop_mode:
                    ptype = 'list_complex'
                else:
                    ptype = 'list_double'
                self.declaration.add((ptype,'P%s' % self.outgoing))
            else:
                coeff = 'COUP'
            for ind in numerator.listindices():
                out.write('  %s[%d] = %s * %s;\n' % (self.outname, 
                                        self.pass_to_HELAS(ind), coeff,
                                        self.write_obj(numerator.get_rep(ind))))
        return out.getvalue()

    def write_obj_Add(self, obj, prefactor=True):
        """Turns addvariable into a string"""
        data = defaultdict(list)
        number = []
        [data[p.prefactor].append(p) if hasattr(p, 'prefactor') else number.append(p)
             for p in obj]
        file_str = StringIO()
        if prefactor and obj.prefactor != 1:
            formatted = self.change_number_format(obj.prefactor)
            if formatted.startswith(('+','-')):
                file_str.write('(%s)' % formatted)
            else:
                file_str.write(formatted)
            file_str.write(' * (')
        else:
            file_str.write('(')
        first=True
        for value, obj_list in data.items():
            add= ' + '
            if value not in  [-1,1]:
                nb_str = self.change_number_format(value)
                if nb_str[0] in ['+','-']:
                    file_str.write(nb_str) # eventually (' '+nb_str)?
                else:
                    file_str.write('+' if first else ' + ')
                    file_str.write(nb_str)
                file_str.write('*(')
            elif value == -1:
                add = ' - ' 
                file_str.write('-' if first else ' - ')
            elif not first:
                file_str.write(' + ')
            else:
                file_str.write('')
            first = False
            file_str.write(add.join([self.write_obj(obj, prefactor=False) 
                                                          for obj in obj_list]))
            if value not in [1,-1]:
                file_str.write(')')
        if number:
            total = sum(number)
            file_str.write('+ %s' % self.change_number_format(total))
        file_str.write(')')
        ###print(file_str.getvalue()) # FOR DEBUGGING
        return file_str.getvalue()


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
        replace_dict['info_lines'] = export_cpp.get_mg5_info_lines().replace('# ','//')
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

    def get_process_function_definitions(self, write=True):
        """The complete Pythia 8 class definition for the process"""
        replace_dict = super(export_cpp.OneProcessExporterGPU,self).get_process_function_definitions(write=False)
        replace_dict['ncouplings'] = len(self.couplings2order)
        replace_dict['ncouplingstimes2'] = 2 *  replace_dict['ncouplings']
        replace_dict['nparams'] = len(self.params2order)
        replace_dict['nmodels'] = replace_dict['nparams'] + replace_dict['ncouplings']
        replace_dict['coupling_list'] = ' '
        coupling = [''] * len(self.couplings2order)
        params = [''] * len(self.params2order)
        for coup, pos in self.couplings2order.items():
            coupling[pos] = coup
        coup_str = "static cxtype tIPC[%s] = {cxmake(pars->%s)};\n"\
            %(len(self.couplings2order), '),cxmake(pars->'.join(coupling))
        for para, pos in self.params2order.items():
            params[pos] = para            
        param_str = "static fptype tIPD[%s] = {(fptype)pars->%s};\n"\
            %(len(self.params2order), ',(fptype)pars->'.join(params))            
        replace_dict['assign_coupling'] = coup_str + param_str
        replace_dict['all_helicities'] = self.get_helicity_matrix(self.matrix_elements[0])
        replace_dict['all_helicities'] = replace_dict['all_helicities'] .replace("helicities", "tHel")
        file = self.read_template_file(self.process_definition_template) %\
               replace_dict
        return file
