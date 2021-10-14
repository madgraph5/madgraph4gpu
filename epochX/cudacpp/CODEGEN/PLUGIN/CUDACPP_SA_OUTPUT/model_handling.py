import os
pjoin = os.path.join

from collections import defaultdict
from fractions import Fraction
from six import StringIO

# AV - use templates for source code, scripts and Makefiles from PLUGINDIR instead of MG5DIR
###from madgraph import MG5DIR
PLUGINDIR = os.path.dirname( __file__ )

# AV - create a plugin-specific logger
import logging
logger = logging.getLogger('madgraph.PLUGIN.CUDACPP_SA_OUTPUT.model_handling')

#------------------------------------------------------------------------------------

# AV - modify export_cpp.get_mg5_info_lines (replace '# ' by '//')
import madgraph.iolibs.export_cpp as export_cpp

def PLUGIN_get_mg5_info_lines():
    return DEFAULT_get_mg5_info_lines().replace('# ','//')    

DEFAULT_get_mg5_info_lines = export_cpp.get_mg5_info_lines
export_cpp.get_mg5_info_lines = PLUGIN_get_mg5_info_lines

#------------------------------------------------------------------------------------

# AV - modify writers.FileWriter.__init__ (add a debug printout)
import madgraph.iolibs.file_writers as writers

def PLUGIN_FileWriter__init__( self, name, opt = 'w' ):
    print( 'FileWriter %s for %s'%( type(self), name) )
    return DEFAULT_FileWriter__init__( self, name, opt )

DEFAULT_FileWriter__init__ = writers.FileWriter.__init__
writers.FileWriter.__init__ = PLUGIN_FileWriter__init__

#------------------------------------------------------------------------------------

# AV - replace writers.CPPWriter by PLUGIN_FileWriter (remove formatting)
class PLUGIN_FileWriter(writers.FileWriter):
    """Default FileWriter with minimal modifications"""

DEFAULT_CPPWriter = writers.CPPWriter
###writers.CPPWriter = DEFAULT_FileWriter # WITH FORMATTING
writers.CPPWriter = PLUGIN_FileWriter # WITHOUT FORMATTING

#------------------------------------------------------------------------------------

import aloha
import aloha.aloha_writers as aloha_writers

class PLUGIN_ALOHAWriter(aloha_writers.ALOHAWriterForGPU):
    # Class structure information
    #  - object
    #  - WriteALOHA(object) [in aloha/aloha_writers.py]
    #  - ALOHAWriterForCPP(WriteALOHA) [in aloha/aloha_writers.py]
    #  - ALOHAWriterForGPU(ALOHAWriterForCPP) [in aloha/aloha_writers.py]
    #  - PLUGIN_ALOHAWriter(ALOHAWriterForGPU)
    #      This class
    
    # AV - keep defaults from aloha_writers.ALOHAWriterForGPU
    ###extension = '.cu'
    ###prefix ='__device__'
    type2def = {}
    type2def['pointer_vertex'] = '*' # using complex<double>* vertex
    type2def['pointer_coup'] = ''

    # AV - modify C++ code from aloha_writers.ALOHAWriterForGPU
    ###ci_definition = 'cxtype cI = cxtype(0., 1.);\n'
    ci_definition = 'const cxtype cI = cxmake( 0., 1. );\n'
    ###realoperator = '.real()'
    ###imagoperator = '.imag()'
    realoperator = 'cxreal' # NB now a function
    imagoperator = 'cximag' # NB now a function

    # AV - improve formatting
    ###type2def['int'] = 'int '
    type2def['int'] = 'int'
    ###type2def['double'] = 'fptype '
    type2def['double'] = 'fptype'
    ###type2def['complex'] = 'cxtype '
    type2def['complex'] = 'cxtype'

    # AV - modify aloha_writers.ALOHAWriterForCPP method (improve formatting)
    def change_number_format(self, number):
        """Formatting the number"""
        def isinteger(x):
            try:
                return int(x) == x
            except TypeError:
                return False
        if isinteger(number):
	    ###out = '%s.' % (str(int(number)))
            if number == 1:
                out = 'one' # AV
            elif number == -1:
                out = '-one' # AV
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

    # AV - modify aloha_writers.ALOHAWriterForCPP method (improve formatting)
    # [NB: this exists in ALOHAWriterForGPU but essentially falls back to ALOHAWriterForCPP]
    # [NB: no, actually this exists twice(!) in ForGPU and the 2nd version is not trivial! but I keep the ForCPP version]
    # This affects HelAmps_sm.h and HelAmps_sm.cu
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
        comment_inputs = [] # AV
        for format, argname in self.define_argument_list(couplings):
            if format.startswith('list'):
                type = self.type2def[format[5:]]
                list_arg = '[]'
                if not argname.startswith('COUP'): type += '_sv' # AV vectorize
                comment_inputs.append('%s[6]'%argname) # AV (wavefuncsize=6 is hardcoded also in export_cpp...)
            else:
                type = self.type2def[format]
                list_arg = ''
            if argname.startswith('COUP'):
                point = self.type2def['pointer_coup']
                ###args.append('%s%s%s%s'% (type, point, argname, list_arg))
                args.append('%s %s%s%s'% (type, point, argname, list_arg)) # AV
            else:
                ###args.append('%s%s%s'% (type, argname, list_arg))
                args.append('%s %s%s'% (type, argname, list_arg)) # AV
        if not self.offshell:
            ###output = '%(doublec)s %(pointer_vertex)s vertex' % {
            output = '%(doublec)s_sv%(pointer_vertex)s vertex' % { # AV vectorize
                'doublec':self.type2def['complex'],
                'pointer_vertex': self.type2def['pointer_vertex']}
            comment_output = 'amplitude \'vertex\''
        else:
            ###output = '%(doublec)s %(spin)s%(id)d[]' % {
            output = '%(doublec)s_sv %(spin)s%(id)d[]' % { # AV vectorize
                     'doublec': self.type2def['complex'],
                     'spin': self.particles[self.outgoing -1],
                     'id': self.outgoing}
            self.declaration.add(('list_complex', output))
            comment_output = 'wavefunction \'%s%d[6]\'' % ( self.particles[self.outgoing -1], self.outgoing ) # AV (wavefuncsize=6)
        ###out.write('%(prefix)s void %(name)s(%(args)s,%(output)s)' % \
        ###          {'prefix': self.prefix,
        ###              'output':output, 'name': name, 'args': ', '.join(args)})
        #out.write('  %(prefix)s void %(name)s( const %(args)s, %(output)s )' % \
        #          {'prefix': self.prefix,
        #              'output':output, 'name': name, 'args': ', const '.join(args)}) # AV - add const
        comment = '// Compute the output %s from the input wavefunctions %s' % ( comment_output, ', '.join(comment_inputs) ) # AV
        indent = ' ' * len( '  void %s( ' % name )
        out.write('  %(comment)s\n  %(prefix)s\n  void %(name)s( const %(args)s,\n%(indent)s%(output)s )%(suffix)s' %
                  {'comment': comment, # AV - add comment
                   'prefix': self.prefix + ( ' INLINE' if 'is_h' in mode else '' ), # AV - add INLINE
                   'suffix': ( ' ALWAYS_INLINE' if 'is_h' in mode else '' ), # AV - add ALWAYS_INLINE
                   'indent':indent, 'output':output, 'name': name,
                   'args': (',\n' + indent + 'const ').join(args)}) # AV - add const, add indent
        if 'is_h' in mode:
            out.write(';\n')
            out.write('\n  //--------------------------------------------------------------------------\n') # AV add footer
        else:
            ###out.write('\n{\n')
            out.write('\n  {\n') # AV
        return out.getvalue() 

    # AV - modify aloha_writers.ALOHAWriterForCPP method (improve formatting)
    # This affects HelAmps_sm.cu
    def get_foot_txt(self):
        """Prototype for language specific footer"""
        ###return '}\n'
        return '  }\n\n  //--------------------------------------------------------------------------' # AV

    # AV - modify aloha_writers.ALOHAWriterForCPP method (improve formatting)
    # This affects HelAmps_sm.cu
    def get_declaration_txt(self, add_i=True):
        """ Prototype for how to write the declaration of variable
            Include the symmetry line (entry FFV_2)
        """        
        out = StringIO()
        out.write('    mgDebug( 0, __FUNCTION__ );\n') # AV - NO! move to get_declaration.txt
        argument_var = [name for type,name in self.call_arg]
        # define the complex number CI = 0+1j
        if add_i:
            ###out.write(self.ci_definition)
            out.write('    ' + self.ci_definition) # AV
        for type, name in self.declaration.tolist():
            if type.startswith('list'):
                type = type[5:]
                if name.startswith('P'):
                    size = 4
                elif not 'tmp' in name:
                    continue # should be defined in the header
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
                ###out.write('    %s %s[%s];\n' % (self.type2def[type], name, size))
                out.write('    %s_sv %s[%s];\n' % (self.type2def[type], name, size)) # AV vectorize
            elif (type, name) not in self.call_arg:
                ###out.write('    %s %s;\n' % (self.type2def[type], name))
                out.write('    %s_sv %s;\n' % (self.type2def[type], name)) # AV vectorize
        ###out.write('    // END DECLARATION\n') # FOR DEBUGGING
        return out.getvalue()

    # AV - modify aloha_writers.ALOHAWriterForCPP method (improve formatting)
    # This affects 'V1[0] = ' in HelAmps_sm.cu
    def get_momenta_txt(self):
        """Define the Header of the C++ file. This include
            - momentum conservation
            - definition of the impulsion"""
        out = StringIO()
        # Define all the required momenta
        p = [] # a list for keeping track how to write the momentum
        signs = self.get_momentum_conservation_sign()
        for i,type in enumerate(self.particles):
            if self.declaration.is_used('OM%s' % (i+1)):
                ###out.write("    OM{0} = {1};\n    if (M{0} != {1})\n OM{0}={2}/(M{0}*M{0});\n".format(
                ###out.write("    OM{0} = {1};\n    if ( M{0} != {1} ) OM{0} = {2} / (M{0}*M{0});\n".format( # AV older
                out.write("    OM{0} = ( M{0} != {1} ? {2} / ( M{0} * M{0} ) : {1} );\n".format( # AV use ternary in OM3
                    ###i+1, self.change_number_format(0), self.change_number_format(1)))
                    i+1, '0.', '1.')) # AV force scalar "1." instead of vector "one"
            if i+1 == self.outgoing:
                out_type = type
                out_size = self.type_to_size[type] 
                continue
            elif self.offshell:
                if len(p) != 0 : p.append(' ') # AV
                ###p.append('{0}{1}{2}[%(i)s]'.format(signs[i],type,i+1,type))
                p.append('{0} {1}{2}[%(i)s]'.format(signs[i],type,i+1,type)) # AV
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
                ###out.write('    %s%s[%s] = %s;\n' % (type,self.outgoing, i, ''.join(p) % dict_energy))
                out.write('    %s%s[%s] = %s;\n' % (type,self.outgoing, i, ''.join(p) % dict_energy)) # AV
            if self.declaration.is_used('P%s' % self.outgoing):
                self.get_one_momenta_def(self.outgoing, out)
        # Returning result
        return out.getvalue()

    # AV - modify aloha_writers.ALOHAWriterForCPP method (improve formatting)
    # This affects 'P1[0] = ' in HelAmps_sm.cu
    def get_one_momenta_def(self, i, strfile):
        type = self.particles[i-1]
        if aloha.loop_mode:
            ###template ='P%(i)d[%(j)d] = %(sign)s%(type)s%(i)d[%(nb)d];\n'
            template ='    P%(i)d[%(j)d] = %(sign)s %(type)s%(i)d[%(nb)d];\n' # AV
            template ='    P%(i)d[%(j)d] = %(sign)s %(type)s%(i)d[%(nb)d];\n' # AV
        else:
            ###template ='P%(i)d[%(j)d] = %(sign)s%(type)s%(i)d[%(nb2)d]%(operator)s;\n'
            ###template ='    P%(i)d[%(j)d] = %(sign)s %(type)s%(i)d[%(nb2)d]%(operator)s;\n' # AV older
            template ='    P%(i)d[%(j)d] = %(sign)s %(operator)s( %(type)s%(i)d[%(nb2)d] );\n' # AV cxreal/cximag
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
	    ###strfile.write(template % {'j':j,'type': type, 'i': i, 
            ###            'nb': nb, 'nb2': nb2, 'operator':operator,
            ###            'sign': self.get_P_sign(i)})
            sign = self.get_P_sign(i) if self.get_P_sign(i) else '+' # AV
            strfile.write(template % {'j':j,'type': type, 'i': i, 
                        'nb': nb, 'nb2': nb2, 'operator':operator,
                        'sign': sign}) # AV

    # AV - modify aloha_writers.ALOHAWriterForCPP method (improve formatting)
    # This is called once per FFV function, i.e. once per WriteALOHA instance?
    # It is called by WriteALOHA.write, after get_header_txt, get_declaration_txt, get_momenta_txt, before get_foot_txt
    # This affects 'denom = COUP' in HelAmps_sm.cu
    # This affects 'V1[2] = ' and 'F1[2] = ' in HelAmps_sm.cu
    # This affects 'TMP0 = ' in HelAmps_sm.cu
    # This affects '(*vertex) = ' in HelAmps_sm.cu
    def define_expression(self):
        """Write the helicity amplitude in C++ format"""
        out = StringIO()
        ###out.write('    mgDebug( 0, __FUNCTION__ );\n') # AV - NO! move to get_declaration.txt
        if self.routine.contracted:
            keys = sorted(self.routine.contracted.keys())
            for name in keys:
                obj = self.routine.contracted[name]
                # This affects 'TMP0 = ' in HelAmps_sm.cu
                ###out.write(' %s = %s;\n' % (name, self.write_obj(obj)))
                out.write('    %s = %s;\n' % (name, self.write_obj(obj))) # AV
                self.declaration.add(('complex', name))
        for name, (fct, objs) in self.routine.fct.items():
            format = ' %s = %s;\n' % (name, self.get_fct_format(fct))
            out.write(format % ','.join([self.write_obj(obj) for obj in objs])) # AV not used in eemumu?
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
                # This affects '(*vertex) = ' in HelAmps_sm.cu
                ###out.write(' %(pre_vertex)svertex%(post_vertex)s = %(pre_coup)sCOUP%(post_coup)s*%(num)s;\n' % mydict)
                out.write('    %(pre_vertex)svertex%(post_vertex)s = %(pre_coup)sCOUP%(post_coup)s * %(num)s;\n' % mydict) # AV
            else:
                mydict= {}
                if self.type2def['pointer_vertex'] in ['*']:
                    mydict['pre_vertex'] = '(*'
                    mydict['post_vertex'] = ')'
                else:
                    mydict['pre_vertex'] = ''
                    mydict['post_vertex'] = ''                 
                mydict['data'] = self.write_obj(numerator.get_rep([0]))
                # This affects '(*vertex) = ' in HelAmps_sm.cu
                ###out.write(' %(pre_vertex)svertex%(post_vertex)s = %(data)s;\n' % mydict)
                out.write('    %(pre_vertex)svertex%(post_vertex)s = %(data)s;\n' % mydict) # AV
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
                        # This affects 'denom = COUP' in HelAmps_sm.cu
                        ###out.write('    denom = %(pre_coup)s%(coup)s%(post_coup)s/(%(denom)s)\n' % mydict) 
                        out.write('    denom = %(pre_coup)s%(coup)s%(post_coup)s / (%(denom)s)\n' % mydict) # AV
                    else:
                        # This affects 'denom = COUP' in HelAmps_sm.cu
                        ###out.write('    denom = %(pre_coup)s%(coup)s%(post_coup)s/((P%(i)s[0]*P%(i)s[0])-(P%(i)s[1]*P%(i)s[1])-(P%(i)s[2]*P%(i)s[2])-(P%(i)s[3]*P%(i)s[3]) - M%(i)s * (M%(i)s -cI* W%(i)s));\n' % mydict)
                        out.write('    denom = %(pre_coup)s%(coup)s%(post_coup)s / ((P%(i)s[0]*P%(i)s[0]) - (P%(i)s[1]*P%(i)s[1]) - (P%(i)s[2]*P%(i)s[2]) - (P%(i)s[3]*P%(i)s[3]) - M%(i)s*(M%(i)s-cI*W%(i)s));\n' % mydict) # AV
                else:
                    if self.routine.denominator:
                        raise Exception('modify denominator are not compatible with complex mass scheme')                
                    # This affects 'denom = COUP' in HelAmps_sm.cu
                    ###out.write('    denom = %(pre_coup)s%(coup)s%(post_coup)s/((P%(i)s[0]*P%(i)s[0])-(P%(i)s[1]*P%(i)s[1])-(P%(i)s[2]*P%(i)s[2])-(P%(i)s[3]*P%(i)s[3]) - (M%(i)s*M%(i)s));\n' % mydict)
                    out.write('    denom = %(pre_coup)s%(coup)s%(post_coup)s / ((P%(i)s[0]*P%(i)s[0])-(P%(i)s[1]*P%(i)s[1])-(P%(i)s[2]*P%(i)s[2])-(P%(i)s[3]*P%(i)s[3]) - (M%(i)s*M%(i)s));\n' % mydict) # AV
                self.declaration.add(('complex','denom'))
                if aloha.loop_mode:
                    ptype = 'list_complex'
                else:
                    ptype = 'list_double'
                self.declaration.add((ptype,'P%s' % self.outgoing))
            else:
                coeff = 'COUP'
            for ind in numerator.listindices():
                # This affects 'V1[2] = ' and 'F1[2] = ' in HelAmps_sm.cu
                ###out.write('    %s[%d]= %s*%s;\n' % (self.outname,
                out.write('    %s[%d] = %s * %s;\n' % (self.outname, # AV
                                        self.pass_to_HELAS(ind), coeff,
                                        self.write_obj(numerator.get_rep(ind))))
        out.write('    mgDebug( 1, __FUNCTION__ );\n') # AV
        out.write('    return;\n') # AV
        return out.getvalue()

    # AV - modify aloha_writers.WriteALOHA method (improve formatting)
    # This affects 'V1[2] = ' and 'F1[2] = ' in HelAmps_sm.cu
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
            ###file_str.write('*(')
            file_str.write(' * (') # AV
        else:
            file_str.write('(')
        first=True
        for value, obj_list in data.items():
            add= ' + '
            if value not in  [-1,1]:
                nb_str = self.change_number_format(value)
                if nb_str[0] in ['+','-']:
                    file_str.write(nb_str) # AV - eventually (' '+nb_str)?
                else:
                    ###file_str.write('+')
                    file_str.write('+' if first else ' + ') # AV
                    file_str.write(nb_str)
                file_str.write('*(')
            elif value == -1:
                ###add = '-'
                ###file_str.write('-')
                add = ' - ' # AV
                file_str.write('-' if first else ' - ') # AV
            elif not first:
                ###file_str.write('+')
                file_str.write(' + ') # AV
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
        ###print(file_str.getvalue()) # AV - FOR DEBUGGING
        return file_str.getvalue()

#------------------------------------------------------------------------------------

class PLUGIN_UFOModelConverter(export_cpp.UFOModelConverterGPU):
    # Class structure information
    #  - object
    #  - UFOModelConverterCPP(object) [in madgraph/iolibs/export_cpp.py]
    #  - UFOModelConverterGPU(UFOModelConverterCPP) [in madgraph/iolibs/export_cpp.py]
    #  - PLUGIN_UFOModelConverter(UFOModelConverterGPU)
    #      This class

    # AV - keep defaults from export_cpp.UFOModelConverterCPP
    ###include_dir = '.'
    ###c_file_dir = '.'
    ###param_template_h = 'cpp_model_parameters_h.inc'
    ###param_template_cc = 'cpp_model_parameters_cc.inc'

    # AV - change defaults from export_cpp.UFOModelConverterCPP
    # (custom tag to appear in 'This file has been automatically generated for')
    output_name = 'CUDA/C++ standalone'

    # AV - change defaults from export_cpp.UFOModelConverterGPU
    ###cc_ext = 'cu' # create HelAmps_sm.cu
    cc_ext = 'cc' # create HelAmps_sm.cc

    # AV - keep defaults from export_cpp.UFOModelConverterGPU
    ###cc_ext = 'cu'
    ###aloha_template_h = pjoin('gpu','cpp_hel_amps_h.inc')
    ###aloha_template_cc = pjoin('gpu','cpp_hel_amps_cc.inc')
    ###helas_h = pjoin('gpu', 'helas.h')
    ###helas_cc = pjoin('gpu', 'helas.cu')

    # AV - use a custom ALOHAWriter (NB: this is an argument to WriterFactory.__new__, either a string or a class!)
    ###aloha_writer = 'cudac' # WriterFactory will use ALOHAWriterForGPU
    aloha_writer = PLUGIN_ALOHAWriter # WriterFactory will use ALOHAWriterForGPU

    # AV - use template files from PLUGINDIR instead of MG5DIR
    def read_aloha_template_files(self, ext):
        """Read all ALOHA template files with extension ext, strip them of
        compiler options and namespace options, and return in a list"""
        ###path = pjoin(MG5DIR, 'aloha','template_files')
        path = pjoin(PLUGINDIR, 'aloha', 'template_files')
        out = []
        if ext == 'h':
            out.append(open(pjoin(path, self.helas_h)).read())
        else:
            out.append(open(pjoin(path, self.helas_cc)).read())
        return out

    # AV - use the plugin's PLUGIN_OneProcessExporter template_path and __template_path (for aloha_template_h/cc)
    @classmethod
    def read_template_file(cls, filename, classpath=False):
        """Open a template file and return the contents."""
        ###return OneProcessExporterCPP.read_template_file(filename, classpath)
        return PLUGIN_OneProcessExporter.read_template_file(filename, classpath)

    # AV - overload export_cpp.UFOModelConverterCPP method (improve formatting)
    def write_parameters(self, params):
        res = super().write_parameters(params)
        if res == '' : res = '  // (none)'
        else : res = '  ' + res # add leading '  ' after the '// Model' line
        res = res.replace('\n','\n  ')
        res = res.replace(',',', ')
        return res

    # AV - overload export_cpp.UFOModelConverterCPP method (improve formatting)
    def write_set_parameters(self, params):
        res = super().write_set_parameters(params)
        if res == '' : res = '// (none)'
        res = res.replace('\n','\n  ')
        return res

    # AV - overload export_cpp.UFOModelConverterCPP method (improve formatting)
    def write_print_parameters(self, params):
        res = super().write_print_parameters(params)
        if res == '' : res = '// (none)'
        res = res.replace('\n','\n  ')
        return res

    # AV - overload export_cpp.UFOModelConverterCPP method (improve formatting)
    def generate_parameters_class_files(self):
        file_h, file_cc = super().generate_parameters_class_files()
        file_h = file_h[:-1] # remove extra trailing '\n'
        file_cc = file_cc[:-1] # remove extra trailing '\n'
        # [NB: there is a minor bug in export_cpp.UFOModelConverterCPP.generate_parameters_class_files
        # ['independent_couplings' contains dependent parameters, 'dependent parameters' contains independent_couplings]
        # [This only affects the order in which they are printed out - which is now reversed in the templates]
        return file_h, file_cc

    # AV - replace export_cpp.UFOModelConverterCPP method (add explicit std namespace)
    def write_print_parameters(self, params):
        """Write out the lines of independent parameters"""
        # For each parameter, write name = expr;
        res_strings = []
        for param in params:
            ###res_strings.append("cout << setw(20) << \"%s \" << \"= \" << setiosflags(ios::scientific) << setw(10) << %s << endl;" % (param.name, param.name))
            res_strings.append("std::cout << std::setw(20) << \"%s \" << \"= \" << std::setiosflags(std::ios::scientific) << std::setw(10) << %s << std::endl;" % (param.name, param.name)) # AV
        ##return "\n".join(res_strings)
        return "\n  ".join(res_strings) # AV (why was this not necessary before?)

    # AV - replace export_cpp.UFOModelConverterCPP method (add debug printouts)
    # (This is where the loop over FFV functions takes place - I had a hard time to understand it)
    # (Note also that write_combined_cc seems to never be called for our eemumu and ggttgg examples)
    # The calling sequence is the following (understood via MG5_debug after forcing an error by renaming 'write')
    # - madgraph_interface.py 8369 in finalize => self._curr_exporter.convert_model(self._curr_model
    # - output.py 127 in convert_model => super().convert_model(model, wanted_lorentz, wanted_coupling)
    # - export_cpp.py 2503 in convert_model => model_builder.write_files()
    # - export_cpp.py 128 in write_files => self.write_aloha_routines()
    # - export_cpp.py 392 in write_aloha_routines => h_rout, cc_rout = abstracthelas.write(output_dir=None,
    # - create_aloha.py 97 in write => text = writer.write(mode=mode, **opt)
    #   [this is PLUGIN_ALOHAWriter.write which defaults to ALOHAWriterForCPP.write]
    #   [therein, cc_text comes from WriteALOHA.write, while h_text comes from get_h_text]
    def write_aloha_routines(self):
        """Generate the hel_amps_model.h and hel_amps_model.cc files, which
        have the complete set of generalized Helas routines for the model"""        
        import aloha.create_aloha as create_aloha
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
        # Read in the template .h and .cc files, stripped of compiler commands and namespaces
        template_h_files = self.read_aloha_template_files(ext = 'h')
        template_cc_files = self.read_aloha_template_files(ext = 'cc')
        aloha_model = create_aloha.AbstractALOHAModel(self.model.get('name'), explicit_combine=True)
        aloha_model.add_Lorentz_object(self.model.get('lorentz'))
        if self.wanted_lorentz:
            aloha_model.compute_subset(self.wanted_lorentz)
        else:
            aloha_model.compute_all(save=False, custom_propa=True)
        for abstracthelas in dict(aloha_model).values():
            print(type(abstracthelas), abstracthelas.name) # AV this is the loop on FFV functions
            h_rout, cc_rout = abstracthelas.write(output_dir=None, language=self.aloha_writer, mode='no_include')
            template_h_files.append(h_rout)
            template_cc_files.append(cc_rout)
        replace_dict['function_declarations'] = '\n'.join(template_h_files)
        replace_dict['function_definitions'] = '\n'.join(template_cc_files)
        file_h = self.read_template_file(self.aloha_template_h) % replace_dict
        file_cc = self.read_template_file(self.aloha_template_cc) % replace_dict
        # Write the files
        writers.CPPWriter(model_h_file).writelines(file_h)
        writers.CPPWriter(model_cc_file).writelines(file_cc)
        logger.info("Created files %s and %s in directory" \
                    % (os.path.split(model_h_file)[-1],
                       os.path.split(model_cc_file)[-1]))
        logger.info("%s and %s" % \
                    (os.path.split(model_h_file)[0],
                     os.path.split(model_cc_file)[0]))

#------------------------------------------------------------------------------------

import madgraph.iolibs.files as files
import madgraph.various.misc as misc

class PLUGIN_OneProcessExporter(export_cpp.OneProcessExporterGPU):
    # Class structure information
    #  - object
    #  - OneProcessExporterCPP(object) [in madgraph/iolibs/export_cpp.py]
    #  - OneProcessExporterGPU(OneProcessExporterCPP) [in madgraph/iolibs/export_cpp.py]
    #  - PLUGIN_OneProcessExporter(OneProcessExporterGPU)
    #      This class
    
    # AV - change defaults from export_cpp.OneProcessExporterGPU
    # [NB process_class = "CPPProcess" is set in OneProcessExporterCPP.__init__]
    # [NB process_class = "gCPPProcess" is set in OneProcessExporterGPU.__init__]
    ###cc_ext = 'cu' # create gCPPProcess.cu (and symlink it as CPPProcess.cc)
    cc_ext = 'cc' # create CPPProcess.cc (and symlink it as gCPPProcess.cu)

    # AV - keep defaults from export_cpp.OneProcessExporterGPU
    ###process_dir = '.'
    ###include_dir = '.'
    ###process_template_h = 'gpu/process_h.inc'
    ###process_template_cc = 'gpu/process_cc.inc'
    ###process_class_template = 'gpu/process_class.inc'
    ###process_definition_template = 'gpu/process_function_definitions.inc'
    ###process_wavefunction_template = 'cpp_process_wavefunctions.inc'
    ###process_sigmaKin_function_template = 'gpu/process_sigmaKin_function.inc'
    ###single_process_template = 'gpu/process_matrix.inc'

    # AV - use template files from PLUGINDIR instead of MG5DIR
    ###template_path = os.path.join(_file_path, 'iolibs', 'template_files')
    ###__template_path = os.path.join(_file_path, 'iolibs', 'template_files') 
    template_path = os.path.join( PLUGINDIR, 'madgraph', 'iolibs', 'template_files' )
    __template_path = os.path.join( PLUGINDIR, 'madgraph', 'iolibs', 'template_files' )

    # AV - overload export_cpp.OneProcessExporterGPU constructor (rename gCPPProcess to CPPProcess)
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.process_class = "CPPProcess"

    # AV - modify export_cpp.OneProcessExporterGPU method (fix gCPPProcess.cu)
    def get_process_function_definitions(self, write=True):
        """The complete class definition for the process"""
        replace_dict = super(export_cpp.OneProcessExporterGPU,self).get_process_function_definitions(write=False)
        replace_dict['ncouplings'] = len(self.couplings2order)
        replace_dict['ncouplingstimes2'] = 2 * replace_dict['ncouplings']
        replace_dict['nparams'] = len(self.params2order)
        replace_dict['nmodels'] = replace_dict['nparams'] + replace_dict['ncouplings']
        replace_dict['coupling_list'] = ' '
        coupling = [''] * len(self.couplings2order)
        params = [''] * len(self.params2order)
        for coup, pos in self.couplings2order.items():
            coupling[pos] = coup
        ###coup_str = "static cxtype tIPC[%s] = {pars->%s};\n"\
        ###    %(len(self.couplings2order), ',pars->'.join(coupling))
        coup_str = "static cxtype tIPC[%s] = { cxmake(pars->%s) };\n"\
            %(len(self.couplings2order), '), cxmake(pars->'.join(coupling)) # AV
        for para, pos in self.params2order.items():
            params[pos] = para
        ###param_str = "static double tIPD[%s] = {pars->%s};\n"\
        ###    %(len(self.params2order), ',pars->'.join(params))
        param_str = "    static fptype tIPD[%s] = { (fptype)pars->%s };"\
            %(len(self.params2order), ', (fptype)pars->'.join(params)) # AV
        replace_dict['assign_coupling'] = coup_str + param_str
        replace_dict['all_helicities'] = self.get_helicity_matrix(self.matrix_elements[0])
        replace_dict['all_helicities'] = replace_dict['all_helicities'] .replace("helicities", "tHel")
        file = self.read_template_file(self.process_definition_template) % replace_dict
        return file

    # AV - modify export_cpp.OneProcessExporterGPU method (fix gCPPProcess.cu)
    def get_all_sigmaKin_lines(self, color_amplitudes, class_name):
        """Get sigmaKin_process for all subprocesses for gCPPProcess.cu"""
        ret_lines = []
        if self.single_helicities:
            ###ret_lines.append( "__device__ void calculate_wavefunctions(int ihel, const fptype* allmomenta,fptype &meHelSum \n#ifndef __CUDACC__\n                                , const int ievt\n#endif\n                                )\n{" )
            ###ret_lines.append(" using namespace MG5_%s;" % self.model_name)
            ###ret_lines.append("mgDebug( 0, __FUNCTION__ );")
            ###ret_lines.append("cxtype amp[1]; // was %i" % len(self.matrix_elements[0].get_all_amplitudes()))
            ###ret_lines.append("const int ncolor =  %i;" % len(color_amplitudes[0]))
            ###ret_lines.append("cxtype jamp[ncolor];")
            ###ret_lines.append("// Calculate wavefunctions for all processes")
            ###ret_lines.append("using namespace MG5_%s;" % self.model_name)
            ret_lines.append('  __device__ void calculate_wavefunctions( int ihel,')
            indent = ' ' * ( ret_lines[-1].find('(') + 2 )
            ret_lines.append(indent+'const fptype* allmomenta,')
            ret_lines.append(indent+'fptype& meHelSum')
            ret_lines.append('#ifndef __CUDACC__')
            ret_lines.append(indent+', const int ievt')
            ret_lines.append('#endif')
            ret_lines.append(indent+')')
            ret_lines.append('  {')
            ret_lines.append('    using namespace MG5_%s;' % self.model_name)
            ret_lines.append('    mgDebug( 0, __FUNCTION__ );')
            ret_lines.append('    cxtype amp[1]; // was %i' % len(self.matrix_elements[0].get_all_amplitudes()))
            ret_lines.append('    const int ncolor = %i;' % len(color_amplitudes[0]))
            ret_lines.append('    cxtype jamp[ncolor];\n')
            ret_lines.append('    // Calculate wavefunctions for all processes')
            helas_calls = self.helas_call_writer.get_matrix_element_calls(\
                                                    self.matrix_elements[0],
                                                    color_amplitudes[0]
                                                    )
            logger.debug("only one Matrix-element supported?")
            self.couplings2order = self.helas_call_writer.couplings2order
            self.params2order = self.helas_call_writer.params2order
            nwavefuncs = self.matrix_elements[0].get_number_of_wavefunctions()
            ###ret_lines.append("cxtype w[nwf][nw6];")
            ret_lines.append('    cxtype w[nwf][nw6];')
            ret_lines += helas_calls
        else:
            ret_lines.extend([self.get_sigmaKin_single_process(i, me) \
                                  for i, me in enumerate(self.matrix_elements)])
        ###to_add = [] # AV - what is this for? comment it out
        ###to_add.extend([self.get_matrix_single_process(i, me,
        ###                                                 color_amplitudes[i],
        ###                                                 class_name) \
        ###                        for i, me in enumerate(self.matrix_elements)])
        ret_lines.extend([self.get_matrix_single_process(i, me,
                                                         color_amplitudes[i],
                                                         class_name) \
                                for i, me in enumerate(self.matrix_elements)])
        return "\n".join(ret_lines)

    # AV - modify export_cpp.OneProcessExporterGPU method (replace '# Process' by '// Process')
    def get_process_info_lines(self, matrix_element):
        """Return info lines describing the processes for this matrix element"""
        ###return"\n".join([ "# " + process.nice_string().replace('\n', '\n# * ') \
        ###                 for process in matrix_element.get('processes')])
        return"\n".join([ "// " + process.nice_string().replace('\n', '\n// * ') \
                         for process in matrix_element.get('processes')])

    # AV - replace the export_cpp.OneProcessExporterGPU method (invert .cc/.cu, add debug printouts)
    def generate_process_files(self):
        """Generate mgOnGpuConfig.h, CPPProcess.cc, check_sa.cc, gCPPProcess.h/cu link, gcheck_sa.cu link""" 
        misc.sprint('Entering PLUGIN_OneProcessExporter.generate_process_files')
        super(export_cpp.OneProcessExporterGPU, self).generate_process_files()
        self.edit_check_sa()
        self.edit_mgonGPU()
        # Add symbolic links
        ###files.ln(pjoin(self.path, 'gcheck_sa.cu'), self.path, 'check_sa.cc')
        ###files.ln(pjoin(self.path, 'gCPPProcess.cu'), self.path, 'CPPProcess.cc')
        files.ln(pjoin(self.path, 'check_sa.cc'), self.path, 'gcheck_sa.cu')
        files.ln(pjoin(self.path, 'CPPProcess.cc'), self.path, 'gCPPProcess.cu')

    # AV - replace the export_cpp.OneProcessExporterGPU method (invert .cc/.cu, add debug printouts)
    def edit_check_sa(self):
        """Generate check_sa.cc"""
        misc.sprint('Entering PLUGIN_OneProcessExporter.edit_check_sa')
        template = open(pjoin(self.template_path,'gpu','check_sa.cc'),'r').read()
        replace_dict = {}
        replace_dict['nexternal'], _ = self.matrix_elements[0].get_nexternal_ninitial()
        replace_dict['model'] = self.model_name
        replace_dict['numproc'] = len(self.matrix_elements)
        ff = open(pjoin(self.path, 'check_sa.cc'),'w')
        ff.write(template)
        ff.close()

    # AV - add debug printouts over the export_cpp.OneProcessExporterGPU method
    def edit_mgonGPU(self):
        """Generate mgOnGpuConfig.h"""
        misc.sprint('Entering PLUGIN_OneProcessExporter.edit_mgonGPU')
        ###misc.sprint('  template_path=%s'%self.template_path) # look for gpu/mgOnGpuConfig.h here
        return super().edit_mgonGPU()

    # AV - overload the export_cpp.OneProcessExporterGPU method (add debug printout and truncate last \n)
    # [*NB export_cpp.UFOModelConverterGPU.write_process_h_file is not called!*]
    def write_process_h_file(self, writer):
        """Generate final gCPPProcess.h"""
        misc.sprint('Entering PLUGIN_OneProcessExporter.write_process_h_file')
        out = super().write_process_h_file(writer)
        writer.seek(-1, os.SEEK_CUR)
        writer.truncate()
        return out

    # AV - replace the export_cpp.OneProcessExporterGPU method (replace HelAmps.cu by HelAmps.cc)
    def super_write_process_cc_file(self, writer):
        """Write the class member definition (.cc) file for the process described by matrix_element"""
        replace_dict = super(export_cpp.OneProcessExporterGPU, self).write_process_cc_file(False)
        ###replace_dict['hel_amps_def'] = "\n#include \"../../src/HelAmps_%s.cu\"" % self.model_name
        replace_dict['hel_amps_def'] = "\n#include \"../../src/HelAmps_%s.cc\"" % self.model_name # AV
        if writer:
            file = self.read_template_file(self.process_template_cc) % replace_dict
            # Write the file
            writer.writelines(file)
        else:
            return replace_dict

    # AV - overload the export_cpp.OneProcessExporterGPU method (add debug printout and truncate last \n)
    def write_process_cc_file(self, writer):
        """Generate CPPProcess.cc"""
        misc.sprint('Entering PLUGIN_OneProcessExporter.write_process_cc_file')
        ###out = super().write_process_cc_file(writer)
        out = self.super_write_process_cc_file(writer)
        writer.seek(-1, os.SEEK_CUR)
        writer.truncate()
        return out

    # AV - replace the export_cpp.OneProcessExporterGPU method (improve formatting)
    @staticmethod
    def coeff(ff_number, frac, is_imaginary, Nc_power, Nc_value=3):
        """Returns a nicely formatted string for the coefficients in JAMP lines"""
        total_coeff = ff_number * frac * Fraction(Nc_value) ** Nc_power
        if total_coeff == 1:
            if is_imaginary:
                ###return '+cxtype(0,1)*'
                return '+cxtype(0, 1) * ' # AV
            else:
                return '+'
        elif total_coeff == -1:
            if is_imaginary:
                ###return '-cxtype(0,1)*'
                return '-cxtype(0, 1) * ' # AV
            else:
                return '-'
        res_str = '%+i.' % total_coeff.numerator
        if total_coeff.denominator != 1:
            # Check if total_coeff is an integer
            res_str = res_str + '/%i.' % total_coeff.denominator
        if is_imaginary:
            res_str = res_str + '*cxtype(0,1)'    
        return res_str + '*'

    # AV - replace the export_cpp.OneProcessExporterCPP method (fix fptype and improve formatting)
    def get_color_matrix_lines(self, matrix_element):
        """Return the color matrix definition lines for this matrix element. Split rows in chunks of size n."""
        import madgraph.core.color_algebra as color
        if not matrix_element.get('color_matrix'):
            ###return "\n".join(["static const double denom[1] = {1.};", "static const double cf[1][1] = {1.};"])
            return "\n".join(["    static const fptype denom[1] = {1.};", "static const fptype cf[1][1] = {1.};"]) # AV
        else:
            color_denominators = matrix_element.get('color_matrix').\
                                                 get_line_denominators()
            ###denom_string = "static const double denom[ncolor] = {%s};" % ",".join(["%i" % denom for denom in color_denominators])
            denom_string = "    static const fptype denom[ncolor] = {%s};" % ", ".join(["%i" % denom for denom in color_denominators]) # AV
            matrix_strings = []
            my_cs = color.ColorString()
            for index, denominator in enumerate(color_denominators):
                # Then write the numerators for the matrix elements
                num_list = matrix_element.get('color_matrix').get_line_numerators(index, denominator)
                ###matrix_strings.append("{%s}" % ",".join(["%d" % i for i in num_list]))
                matrix_strings.append("{%s}" % ", ".join(["%d" % i for i in num_list])) # AV
            ###matrix_string = "static const double cf[ncolor][ncolor] = {" + ",".join(matrix_strings) + "};"
            matrix_string = "    static const fptype cf[ncolor][ncolor] = {\n      " + ",\n      ".join(matrix_strings) + "};" # AV
            return "\n".join([denom_string, matrix_string])

    # AV - replace the export_cpp.OneProcessExporterGPU method (improve formatting)
    def get_initProc_lines(self, matrix_element, color_amplitudes):
        """Get initProc_lines for function definition for gCPPProcess::initProc"""
        initProc_lines = []
        initProc_lines.append("// Set external particle masses for this matrix element")
        for part in matrix_element.get_external_wavefunctions():
            ###initProc_lines.append("mME.push_back(pars->%s);" % part.get('mass'))
            initProc_lines.append("    mME.push_back( pars->%s );" % part.get('mass')) # AV
        ###for i, colamp in enumerate(color_amplitudes):
        ###    initProc_lines.append("jamp2[%d] = new double[%d];" % (i, len(colamp))) # AV - this was commented out already
        return "\n".join(initProc_lines)

    # AV - replace the export_cpp.OneProcessExporterCPP method (improve formatting)
    def get_helicity_matrix(self, matrix_element):
        """Return the Helicity matrix definition lines for this matrix element"""
        ###helicity_line = "static const int helicities[ncomb][nexternal] = {";
        helicity_line = "    static const int helicities[ncomb][nexternal] = {\n      "; # AV
        helicity_line_list = []
        for helicities in matrix_element.get_helicity_matrix(allow_reverse=False):
            ###helicity_line_list.append("{"+",".join(['%d'] * len(helicities)) % tuple(helicities) + "}")
            helicity_line_list.append( "{" + ", ".join(['%d'] * len(helicities)) % tuple(helicities) + "}" ) # AV
        ###return helicity_line + ",".join(helicity_line_list) + "};"
        return helicity_line + ",\n      ".join(helicity_line_list) + "};" # AV

    # AV - overload the export_cpp.OneProcessExporterGPU method (just to add some comments...)
    def get_reset_jamp_lines(self, color_amplitudes):
        """Get lines to reset jamps"""
        ret_lines = super().get_reset_jamp_lines(color_amplitudes)
        if ret_lines != "" : ret_lines = '    // Reset jamp (reset color flows)\n' + ret_lines # AV THIS SHOULD NEVER HAPPEN!
        return ret_lines

#------------------------------------------------------------------------------------

import madgraph.iolibs.helas_call_writers as helas_call_writers

# AV - define a custom HelasCallWriter
class PLUGIN_GPUFOHelasCallWriter(helas_call_writers.GPUFOHelasCallWriter):
    """ A Custom HelasCallWriter """
    # Class structure information
    #  - object
    #  - dict(object) [built-in]
    #  - PhysicsObject(dict) [in madgraph/core/base_objects.py]
    #  - HelasCallWriter(base_objects.PhysicsObject) [in madgraph/iolibs/helas_call_writers.py]
    #  - UFOHelasCallWriter(HelasCallWriter) [in madgraph/iolibs/helas_call_writers.py]
    #  - CPPUFOHelasCallWriter(UFOHelasCallWriter) [in madgraph/iolibs/helas_call_writers.py]
    #  - GPUFOHelasCallWriter(CPPUFOHelasCallWriter) [in madgraph/iolibs/helas_call_writers.py]
    #  - PLUGIN_GPUFOHelasCallWriter(GPUFOHelasCallWriter)
    #      This class

    # AV - replace helas_call_writers.GPUFOHelasCallWriter method (improve formatting of gCPPProcess.cu)
    # [GPUFOHelasCallWriter.format_coupling is called by GPUFOHelasCallWriter.get_external_line/generate_helas_call]
    # [GPUFOHelasCallWriter.get_external_line is called by GPUFOHelasCallWriter.get_external]
    # [GPUFOHelasCallWriter.get_external (adding #ifdef CUDA) is called by GPUFOHelasCallWriter.generate_helas_call]
    # [GPUFOHelasCallWriter.generate_helas_call is called by UFOHelasCallWriter.get_wavefunction_call/get_amplitude_call]
    ###findcoupling = re.compile('pars->([-]*[\d\w_]+)\s*,')
    def NOTUSED__format_coupling(self, call):
        """Format the coupling so any minus signs are put in front"""
        import re
        model = self.get('model')
        if not hasattr(self, 'couplings2order'):
            self.couplings2order = {}
            self.params2order = {}
        for coup in re.findall(self.findcoupling, call):
            if coup == 'ZERO':
                call = call.replace('pars->ZERO', '0.')
                continue
            sign = '' 
            if coup.startswith('-'):
                sign = '-'
                coup = coup[1:]
            try:
                param = model.get_parameter(coup)
            except KeyError:
                param = False
            if param:   
                alias = self.params2order
                name = "cIPD"
            else: 
                alias = self.couplings2order
                name = "cIPC"
            if coup not in alias:
                alias[coup] = len(alias)
            if name == "cIPD":
                call = call.replace('pars->%s%s' % (sign, coup), 
                                    '%s%s[%s]' % (sign, name, alias[coup]))
            else:
                call = call.replace('pars->%s%s' % (sign, coup), 
                                    '%scxtype(cIPC[%s],cIPC[%s])' % 
                                    (sign, 2*alias[coup],2*alias[coup]+1))
        return call

    # AV - new method for formatting wavefunction/amplitude calls
    # [It would be too complex to modify them in helas_objects.HelasWavefunction/Amplitude.get_call_key]
    @staticmethod
    def format_call(call):
        return call.replace('(','( ').replace(')',' )').replace(',',', ')

    # AV - replace helas_call_writers.GPUFOHelasCallWriter method (improve formatting)
    def super_get_matrix_element_calls(self, matrix_element, color_amplitudes):
        """Return a list of strings, corresponding to the Helas calls for the matrix element"""
        import madgraph.core.helas_objects as helas_objects
        import madgraph.loop.loop_helas_objects as loop_helas_objects
        assert isinstance(matrix_element, helas_objects.HelasMatrixElement), \
               "%s not valid argument for get_matrix_element_calls" % \
               type(matrix_element)
        # Do not reuse the wavefunctions for loop matrix elements
        if isinstance(matrix_element, loop_helas_objects.LoopHelasMatrixElement):
            return self.get_loop_matrix_element_calls(matrix_element)
        # Restructure data for easier handling
        color = {}
        for njamp, coeff_list in enumerate(color_amplitudes):
            for coeff, namp in coeff_list:
                if namp not in color:
                    color[namp] = {}
                color[namp][njamp] = coeff
        me = matrix_element.get('diagrams')
        matrix_element.reuse_outdated_wavefunctions(me)
        res = []
        ###res.append('for(int i=0;i<%s;i++){jamp[i] = cxtype(0.,0.);}' % len(color_amplitudes))
        res.append('for( int i=0; i<%s; i++ ){ jamp[i] = cxtype( 0., 0. ); } // reset jamp (reset color flows)' % len(color_amplitudes)) # AV
        for diagram in matrix_element.get('diagrams'):
            ###print('DIAGRAM %3d: #wavefunctions=%3d, #diagrams=%3d' %
            ###      (diagram.get('number'), len(diagram.get('wavefunctions')), len(diagram.get('amplitudes')) )) # AV - FOR DEBUGGING
            res.append('\n    // *** DIAGRAM %d OF %d ***' % (diagram.get('number'), len(matrix_element.get('diagrams'))) ) # AV
            res.append('\n    // Wavefunction(s) for diagram number %d' % diagram.get('number')) # AV
            ###res.extend([ self.get_wavefunction_call(wf) for wf in diagram.get('wavefunctions') ])
            res.extend([ self.format_call(self.get_wavefunction_call(wf)) for wf in diagram.get('wavefunctions') ]) # AV
            if len(diagram.get('wavefunctions')) == 0 : res.append('// (none)') # AV
            ###res.append("# Amplitude(s) for diagram number %d" % diagram.get('number'))
            res.append("\n    // Amplitude(s) for diagram number %d" % diagram.get('number'))
            for amplitude in diagram.get('amplitudes'):
                namp = amplitude.get('number')
                amplitude.set('number', 1)
                ###res.append(self.get_amplitude_call(amplitude))
                res.append(self.format_call(self.get_amplitude_call(amplitude))) # AV
                for njamp, coeff in color[namp].items():
                    ###res.append("jamp[%s] += %samp[0];" % (njamp, export_cpp.OneProcessExporterGPU.coeff(*coeff)))
                    res.append("jamp[%s] += %samp[0];" % (njamp, PLUGIN_OneProcessExporter.coeff(*coeff)))
            if len(diagram.get('amplitudes')) == 0 : res.append('// (none)') # AV
        ###res.append('\n    // *** END OF DIAGRAMS ***' ) # AV - no longer needed ('COLOR ALGEBRA BELOW')
        return res

    # AV - overload helas_call_writers.GPUFOHelasCallWriter method (improve formatting)
    def get_matrix_element_calls(self, matrix_element, color_amplitudes):
        """Return a list of strings, corresponding to the Helas calls for the matrix element"""
        ###res = super().get_matrix_element_calls(matrix_element, color_amplitudes)
        res = self.super_get_matrix_element_calls(matrix_element, color_amplitudes)
        for i, item in enumerate(res):
            ###print(item) # FOR DEBUGGING
            if item.startswith('# Amplitude'): item='//'+item[1:] # AV replace '# Amplitude' by '// Amplitude'
            if not item.startswith('\n') and not item.startswith('#'): res[i]='    '+item
        return res

    # AV - replace helas_call_writers.GPUFOHelasCallWriter method (improve formatting)
    # [GPUFOHelasCallWriter.format_coupling is called by GPUFOHelasCallWriter.get_external_line/generate_helas_call]
    # [GPUFOHelasCallWriter.get_external_line is called by GPUFOHelasCallWriter.get_external]
    # [=> GPUFOHelasCallWriter.get_external is called by GPUFOHelasCallWriter.generate_helas_call]
    # [GPUFOHelasCallWriter.generate_helas_call is called by UFOHelasCallWriter.get_wavefunction_call/get_amplitude_call]
    def get_external(self,wf, argument):
        ###text = '\n#ifdef __CUDACC__\n    %s    \n#else\n    %s\n#endif \n'
        text = '#ifdef __CUDACC__\n    %s\n#else\n    %s\n#endif\n' # AV
        line = self.get_external_line(wf, argument)
        split_line = line.split(',')
        split_line = [ str.lstrip(' ').rstrip(' ') for str in split_line] # AV
        # (AV join using ',': no need to add a space as this is done by format_call later on)
        line = ','.join(split_line) # AV (for CUDA)
        ###split_line.insert(-1, ' ievt')
        split_line.insert(-1, 'ievt') # AV (for C++)
        return text % (line, ','.join(split_line))
    
    # AV - replace helas_call_writers.GPUFOHelasCallWriter method (improve formatting)
    # [GPUFOHelasCallWriter.get_external_line is called by GPUFOHelasCallWriter.get_external]
    # [GPUFOHelasCallWriter.get_external (adding #ifdef CUDA) is called by GPUFOHelasCallWriter.generate_helas_call]
    # [GPUFOHelasCallWriter.generate_helas_call is called by UFOHelasCallWriter.get_wavefunction_call/get_amplitude_call]
    def NOTUSED__get_external_line(self, wf, argument):
        call = ''
        call = call + helas_call_writers.HelasCallWriter.mother_dict[\
                argument.get_spin_state_number()].lower() 
        if wf.get('mass').lower() != 'zero' or argument.get('spin') != 2: 
            # Fill out with X up to 6 positions
            call = call + 'x' * (6 - len(call))
            # Specify namespace for Helas calls
            ##call = call + "((double *)(dps + %d * dpt),"
            call = call + "(allmomenta,"
            if argument.get('spin') != 1:
                # For non-scalars, need mass and helicity
                call = call + "pars->%s, cHel[ihel][%d],"
            else:
                call = call + "pars->%s,"
            call = call + "%+d,w[%d], %d);"
            if argument.get('spin') == 1:
                return call % \
                                (wf.get('mass'),
                                 # For boson, need initial/final here
                                 (-1) ** (wf.get('state') == 'initial'),
                                 wf.get('me_id')-1,
                                 wf.get('number_external')-1)
            elif argument.is_boson():
                misc.sprint(call)
                misc.sprint( (wf.get('mass'),
                                 wf.get('number_external')-1,
                                 # For boson, need initial/final here
                                 (-1) ** (wf.get('state') == 'initial'),
                                 wf.get('me_id')-1,
                                 wf.get('number_external')-1))
                return  self.format_coupling(call % \
                                (wf.get('mass'),
                                 wf.get('number_external')-1,
                                 # For boson, need initial/final here
                                 (-1) ** (wf.get('state') == 'initial'),
                                 wf.get('me_id')-1,
                                 wf.get('number_external')-1))
            else:
                return self.format_coupling(call % \
                                (wf.get('mass'),
                                 wf.get('number_external')-1,
                                 # For fermions, need particle/antiparticle
                                 - (-1) ** wf.get_with_flow('is_part'),
                                 wf.get('me_id')-1,
                                 wf.get('number_external')-1))
        else:
            if wf.get('number_external') == 1:
                call += 'pz'
            elif wf.get('number_external') == 2:
                call += 'mz'
            else:
                call += 'xz'
            call = call + 'x' * (6 - len(call))
            # Specify namespace for Helas calls
            ##call = call + "((double *)(dps + %d * dpt),"
            call = call + "(allmomenta, cHel[ihel][%d],%+d,w[%d],%d);"
            return self.format_coupling(call % \
                                (wf.get('number_external')-1,
                                 # For fermions, need particle/antiparticle
                                 - (-1) ** wf.get_with_flow('is_part'),
                                 wf.get('me_id')-1,
                                 wf.get('number_external')-1))

# AV - use the custom HelasCallWriter
DEFAULT_GPUFOHelasCallWriter = helas_call_writers.GPUFOHelasCallWriter
helas_call_writers.GPUFOHelasCallWriter = PLUGIN_GPUFOHelasCallWriter

#------------------------------------------------------------------------------------
