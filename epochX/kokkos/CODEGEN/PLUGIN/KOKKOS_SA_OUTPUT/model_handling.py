import os

# AV - use templates for source code, scripts and Makefiles from PLUGINDIR instead of MG5DIR
###from madgraph import MG5DIR
PLUGINDIR = os.path.dirname( __file__ )

# AV - create a plugin-specific logger
import logging
logger = logging.getLogger('madgraph.PLUGIN.KOKKOS_SA_OUTPUT.model_handling')

#------------------------------------------------------------------------------------

# AV - import the independent 2nd copy of the export_cpp module (as PLUGIN_export_cpp), previously loaded in output.py
###import madgraph.iolibs.export_cpp as export_cpp # 1st copy
######import madgraph.iolibs.export_cpp as PLUGIN_export_cpp # this is not enough to define an independent 2nd copy: id(export_cpp)==id(PLUGIN_export_cpp)
import PLUGIN.KOKKOS_SA_OUTPUT.PLUGIN_export_cpp as PLUGIN_export_cpp # 2nd copy loaded in the plugin's output.py
###print('id(export_cpp)=%s'%id(export_cpp))
###print('id(PLUGIN_export_cpp)=%s'%id(PLUGIN_export_cpp))

#------------------------------------------------------------------------------------

# AV - modify export_cpp.get_mg5_info_lines (replace '# ' by '//')
def PLUGIN_get_mg5_info_lines():
    return DEFAULT_get_mg5_info_lines().replace('# ','//')

DEFAULT_get_mg5_info_lines = PLUGIN_export_cpp.get_mg5_info_lines
PLUGIN_export_cpp.get_mg5_info_lines = PLUGIN_get_mg5_info_lines

#------------------------------------------------------------------------------------

# AV - load an independent 2nd copy of the writers module (as PLUGIN_writers) and use that within the plugin (workaround for #341)
# See https://stackoverflow.com/a/11285504
###import madgraph.iolibs.file_writers as writers # 1st copy
import sys
import importlib.util
SPEC_WRITERS = importlib.util.find_spec('madgraph.iolibs.file_writers')
PLUGIN_writers = importlib.util.module_from_spec(SPEC_WRITERS)
SPEC_WRITERS.loader.exec_module(PLUGIN_writers)
###sys.modules['PLUGIN.CUDACPP_SA_OUTPUT.PLUGIN_writers'] = PLUGIN_writers # would allow 'import PLUGIN.CUDACPP_SA_OUTPUT.PLUGIN_writers' (not needed)
del SPEC_WRITERS

# AV - use the independent 2nd copy of the writers module within the PLUGIN_export_cpp module (workaround for #341)
###DEFAULT_writers = PLUGIN_export_cpp.writers # not needed
PLUGIN_export_cpp.writers = PLUGIN_writers

#------------------------------------------------------------------------------------

# AV - modify writers.FileWriter.__init__ (add a debug printout)
def PLUGIN_FileWriter__init__( self, name, opt = 'w' ):
    print( 'FileWriter %s for %s'%( type(self), name) )
    return DEFAULT_FileWriter__init__( self, name, opt )

DEFAULT_FileWriter__init__ = PLUGIN_writers.FileWriter.__init__
PLUGIN_writers.FileWriter.__init__ = PLUGIN_FileWriter__init__

#------------------------------------------------------------------------------------

# AV - replace writers.CPPWriter by PLUGIN_CPPWriter (remove formatting)
class PLUGIN_CPPWriter(PLUGIN_writers.FileWriter):
    """Custom CPPWriter based on the default FileWriter with minimal modifications"""

DEFAULT_CPPWriter = PLUGIN_writers.CPPWriter
###PLUGIN_writers.CPPWriter = DEFAULT_CPPWriter # WITH FORMATTING
PLUGIN_writers.CPPWriter = PLUGIN_CPPWriter # WITHOUT FORMATTING

#------------------------------------------------------------------------------------

import aloha
import aloha.aloha_writers as aloha_writers

from collections import defaultdict
from fractions import Fraction
from six import StringIO

# AV - define a custom ALOHAWriter
# (NB: enable this via PLUGIN_UFOModelConverter.aloha_writer)
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
    extension = '.cc'
    prefix = 'KOKKOS_INLINE_FUNCTION'
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

    # AV - add vector types
    type2def['double_v'] = 'fptype'
    type2def['complex_v'] = 'cxtype'

     # AV - modify C++ code from aloha_writers.ALOHAWriterForGPU
    # AV new option: declare C++ variable type only when they are defined?
    ###nodeclare = False # old behaviour (separate declaration with no initialization)
    nodeclare = True # new behaviour (delayed declaration with initialisation)

    # AV - modify aloha_writers.ALOHAWriterForCPP method (improve formatting)
    def change_number_format(self, number):
        """Formatting the number"""
        def isinteger(x):
            try:
                return int(x) == x
            except TypeError:
                return False
        if isinteger(number):
            if number == 1: out = 'one' # AV
            elif number == -1: out = '-one' # AV
            elif number == 2: out = 'two' # AV
            elif number == -2: out = '-two' # AV
            else: out = '%s.' % (str(int(number))) # This prints -1 as '-1.'
        elif isinstance(number, complex):
            if number.imag:
                if number.real:
                    out = '( %s + %s * cI )' % (self.change_number_format(number.real), \
                                    self.change_number_format(number.imag))
                else:
                    if number.imag == 1:
                        out = 'cI'
                    elif number.imag == -1:
                        out = '-cI'
                    else:
                        out = '( %s * cI )' % self.change_number_format(number.imag)
            else:
                out = '%s' % (self.change_number_format(number.real))
        else:
            tmp = Fraction(str(number))
            tmp = tmp.limit_denominator(100)
            if not abs(tmp - number) / abs(tmp + number) < 1e-8: out = '%.9f' % (number)
            elif tmp.numerator == 1 and tmp.denominator == 2 : out = 'half' # AV
            elif tmp.numerator == -1 and tmp.denominator == 2 : out = '-half' # AV
            else: out = '%s./%s.' % (tmp.numerator, tmp.denominator)
        return out

    # AV - modify aloha_writers.ALOHAWriterForCPP method (improve formatting)
    # [NB: this exists in ALOHAWriterForGPU but essentially falls back to ALOHAWriterForCPP]
    # [NB: no, actually this exists twice(!) in ForGPU and the 2nd version is not trivial! but I keep the ForCPP version]
    # This affects HelAmps_sm.h and HelAmps_sm.cc
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
                type = self.type2def[format[5:]] # double or complex (instead of list_double or list_complex)
                comment_inputs.append('%s[6]'%argname) # AV (wavefuncsize=6 is hardcoded also in export_cpp...)
                ###if not argname.startswith('COUP'): type = self.type2def[format[5:]+'_v'] # AV vectorize (double_v or complex_v)
                if not argname.startswith('COUP'):
                    type = self.type2def['complex'] # AV from cxtype to fptype
                    argname = argname
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
            ###output = '%(doublec)s%(pointer_vertex)s allvertexes' % {
            ###    'doublec': self.type2def['double'],
            ###    'pointer_vertex': self.type2def['pointer_vertex']}
            output = '%(doublec)s vertex[]' % {
                #'%(doublec)s allvertexes[]' % {
                'doublec': self.type2def['complex']}
            comment_output = 'amplitude \'vertex\''
            template = ''
        else:
            output = '%(doublec)s %(spin)s%(id)d[]' % {
                #'%(doublec)s all%(spin)s%(id)d[]' % {
                     'doublec': self.type2def['complex'],
                     'spin': self.particles[self.outgoing -1],
                     'id': self.outgoing}
            ###self.declaration.add(('list_complex', output)) # AV BUG FIX - THIS IS NOT NEEDED AND IS WRONG (adds name 'cxtype V3[]')
            comment_output = 'wavefunction \'%s%d[6]\'' % ( self.particles[self.outgoing -1], self.outgoing ) # AV (wavefuncsize=6)
            template = ''
        comment = '// Compute the output %s from the input wavefunctions %s' % ( comment_output, ', '.join(comment_inputs) ) # AV
        indent = ' ' * len( '  void %s( ' % name )
        out.write('  %(comment)s\n  %(template)s\n  %(prefix)s\n  void %(name)s( const %(args)s,\n%(indent)s%(output)s )%(suffix)s' %
                  {'comment': comment, # AV - add comment
                   'template': template, # AV - add template
                   'prefix':'KOKKOS_INLINE_FUNCTION',
                   'suffix':'',
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
    # This affects HelAmps_sm.cc
    def get_foot_txt(self):
        """Prototype for language specific footer"""
        ###return '}\n'
        return '  }\n\n  //--------------------------------------------------------------------------' # AV

    # AV - modify aloha_writers.ALOHAWriterForCPP method (improve formatting)
    # This affects HelAmps_sm.cc
    def get_declaration_txt(self, add_i=True):
        """ Prototype for how to write the declaration of variable
            Include the symmetry line (entry FFV_2)
        """        
        out = StringIO()
        #out.write('    mgDebug( 0, __FUNCTION__ );\n') # AV - NO! move to get_declaration.txt
        argument_var = [name for type,name in self.call_arg]
        # define the complex number CI = 0+1j
        if add_i:
            ###out.write(self.ci_definition)
            out.write('    ' + self.ci_definition) # AV
        codedict = {} # AV allow delayed declaration with initialisation
        for type, name in self.declaration.tolist():
            ###print(name) # FOR DEBUGGING
            ###out.write('    %s %s;\n' % ( type, name ) ) # FOR DEBUGGING
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
                fullname = '%s[%s]'%(name, size) # AV
            elif (type, name) not in self.call_arg:
                fullname = name # AV
            else:
                continue # AV no need to declare the variable
            if fullname.startswith('OM') :
                codedict[fullname] = '%s %s' % (self.type2def[type], fullname) # AV UGLY HACK (OM3 is always a scalar)
            else:
                codedict[fullname] = '%s %s' % (self.type2def[type+'_v'], fullname) # AV vectorize, add to codedict
            ###print(fullname, codedict[fullname]) # FOR DEBUGGING
            if self.nodeclare:
                self.declaration.codedict = codedict # AV new behaviour (delayed declaration with initialisation)
            else:
                out.write('    %s;\n' % codedict[fullname] ) # AV old behaviour (separate declaration with no initialization)
        ###out.write('    // END DECLARATION\n') # FOR DEBUGGING
        return out.getvalue()

    # AV - modify aloha_writers.ALOHAWriterForCPP method (improve formatting)
    # This affects 'V1[0] = ' in HelAmps_sm.cc
    def get_momenta_txt(self):
        """Define the Header of the C++ file. This include
            - momentum conservation
            - definition of the impulsion"""
        out = StringIO()
        # Define all the required momenta
        p = [] # a list for keeping track how to write the momentum
        signs = self.get_momentum_conservation_sign()
        for i, type in enumerate(self.particles):
            if self.declaration.is_used( 'OM%s' % (i+1) ):
                declname = 'OM%s' % (i+1)
                if self.nodeclare: declname = 'const ' + self.declaration.codedict[declname]
                out.write('    {3} = ( M{0} != {1} ? {2} / ( M{0} * M{0} ) : {1} );\n'.format( # AV use ternary in OM3
                    i+1, '0.', '1.', declname)) # AV force scalar "1." instead of vector "one", add declaration
            if i+1 == self.outgoing:
                out_type = type
                out_size = self.type_to_size[type] 
                continue
            elif self.offshell:
                if len(p) == 0 :
                    p.append('{0}{1}{2}[%(i)s]'.format(signs[i],type,i+1,type)) # AV for clang-format (ugly!)
                else:
                    p.append(' ')
                    p.append('{0} {1}{2}[%(i)s]'.format(signs[i],type,i+1,type))
            if self.declaration.is_used('P%s' % (i+1)):
                self.get_one_momenta_def(i+1, out)
        # Define the resulting momenta
        if self.offshell:
            energy_pos = out_size -2
            type = self.particles[self.outgoing-1]
            if aloha.loop_mode:
                size_p = 4
            else:
                size_p = 2
            for i in range(size_p):
                dict_energy = {'i':i}
                out.write( '    %s%s[%s] = %s;\n' % ( type, self.outgoing, i, ''.join(p) % dict_energy ) )
            if self.declaration.is_used( 'P%s' % self.outgoing ):
                self.get_one_momenta_def( self.outgoing, out )
        # Returning result
        ###print('."' + out.getvalue() + '"') # AV - FOR DEBUGGING
        return out.getvalue()

    # AV - modify aloha_writers.ALOHAWriterForCPP method (improve formatting, add delayed declaration with initialisation)
    # This affects 'P1[0] = ' in HelAmps_sm.cc
    def get_one_momenta_def(self, i, strfile):
        type = self.particles[i-1]
        if aloha.loop_mode:
            ptype = 'complex_v'
            templateval ='%(sign)s %(type)s%(i)d[%(nb)d]' # AV
        else:
            ptype = 'double_v'
            templateval ='%(sign)s%(operator)s( %(type)s%(i)d[%(nb2)d] )' # AV cxreal/cximag
        if self.nodeclare: strfile.write('    const %s P%d[4] = { ' % ( self.type2def[ptype], i) ) # AV
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
            sign = self.get_P_sign(i) if self.get_P_sign(i) else '+' # AV
            if self.nodeclare: template = templateval + ( ', ' if j<3 else '' ) # AV
            else: template ='    P%(i)d[%(j)d] = ' + templateval + ';\n' # AV
            strfile.write(template % {'j':j,'type': type, 'i': i,
                        'nb': nb, 'nb2': nb2, 'operator':operator,
                        'sign': sign}) # AV
        if self.nodeclare: strfile.write(' };\n') # AV

    # AV - modify aloha_writers.ALOHAWriterForCPP method (improve formatting)
    # This is called once per FFV function, i.e. once per WriteALOHA instance?
    # It is called by WriteALOHA.write, after get_header_txt, get_declaration_txt, get_momenta_txt, before get_foot_txt
    # This affects 'denom = COUP' in HelAmps_sm.cc
    # This affects 'V1[2] = ' and 'F1[2] = ' in HelAmps_sm.cc
    # This affects 'TMP0 = ' in HelAmps_sm.cc
    # This affects '( *vertex ) = ' in HelAmps_sm.cc
    def define_expression(self):
        """Write the helicity amplitude in C++ format"""
        out = StringIO()
        ###out.write('    mgDebug( 0, __FUNCTION__ );\n') # AV - NO! move to get_declaration.txt
        if self.routine.contracted:
            keys = sorted(self.routine.contracted.keys())
            for name in keys:
                obj = self.routine.contracted[name]
                # This affects 'TMP0 = ' in HelAmps_sm.cc
                ###out.write(' %s = %s;\n' % (name, self.write_obj(obj)))
                if self.nodeclare:
                    out.write('    const %s %s = %s;\n' %
                              (self.type2def['complex_v'], name, self.write_obj(obj))) # AV
                else:
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
                mydict = {'num': self.write_obj(numerator.get_rep([0]))} # '...(TMP4)-cI...' comes from here
                for c in ['coup', 'vertex']:
                    if self.type2def['pointer_%s' %c] in ['*']:
                        mydict['pre_%s' %c] = '( *'
                        mydict['post_%s' %c] = ' )'
                    else:
                        mydict['pre_%s' %c] = ''
                        mydict['post_%s'%c] = ''
                # This affects '( *vertex ) = ' in HelAmps_sm.cc
                out.write('    %(pre_vertex)svertex%(post_vertex)s = %(pre_coup)sCOUP%(post_coup)s * %(num)s;\n' % mydict)
            else:
                mydict= {}
                if self.type2def['pointer_vertex'] in ['*']:
                    mydict['pre_vertex'] = '( *'
                    mydict['post_vertex'] = ' )'
                else:
                    mydict['pre_vertex'] = ''
                    mydict['post_vertex'] = ''
                mydict['data'] = self.write_obj(numerator.get_rep([0]))
                # This affects '( *vertex ) = ' in HelAmps_sm.cc
                out.write('    %(pre_vertex)svertex%(post_vertex)s = %(data)s;\n' % mydict)
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
                if self.nodeclare:
                    mydict['declnamedenom'] = 'const %s denom' % self.type2def['complex_v'] # AV
                else:
                    mydict['declnamedenom'] = 'denom' # AV
                    self.declaration.add(('complex','denom'))
                if not aloha.complex_mass:
                    # This affects 'denom = COUP' in HelAmps_sm.cc
                    if self.routine.denominator:
                        if self.routine.denominator == '1':
                            out.write('    %(declnamedenom)s = %(pre_coup)s%(coup)s%(post_coup)s\n' % mydict) # AV
                        else:
                            mydict['denom'] = self.routine.denominator
                            out.write('    %(declnamedenom)s = %(pre_coup)s%(coup)s%(post_coup)s / ( %(denom)s )\n' % mydict) # AV
                    else:
                        out.write('    %(declnamedenom)s = %(pre_coup)s%(coup)s%(post_coup)s / ( ( P%(i)s[0] * P%(i)s[0] ) - ( P%(i)s[1] * P%(i)s[1] ) - ( P%(i)s[2] * P%(i)s[2] ) - ( P%(i)s[3] * P%(i)s[3] ) - M%(i)s * ( M%(i)s - cI * W%(i)s ) );\n' % mydict) # AV
                else:
                    if self.routine.denominator:
                        raise Exception('modify denominator are not compatible with complex mass scheme')
                    # This affects 'denom = COUP' in HelAmps_sm.cc
                    out.write('    %(declnamedenom)s = %(pre_coup)s%(coup)s%(post_coup)s / ( ( P%(i)s[0] * P%(i)s[0] ) - ( P%(i)s[1] *P%(i)s[1] ) - ( P%(i)s[2] * P%(i)s[2] ) - ( P%(i)s[3] * P%(i)s[3] ) - ( M%(i)s * M%(i)s ) );\n' % mydict) # AV
                ###self.declaration.add(('complex','denom')) # AV moved earlier (or simply removed)
                if aloha.loop_mode: ptype = 'list_complex'
                else: ptype = 'list_double'
                self.declaration.add((ptype,'P%s' % self.outgoing))
            else:
                coeff = 'COUP'
            for ind in numerator.listindices():
                # This affects 'V1[2] = ' and 'F1[2] = ' in HelAmps_sm.cc
                ###out.write('    %s[%d]= %s*%s;\n' % (self.outname,
                out.write('    %s[%d] = %s * %s;\n' % (self.outname, # AV
                                        self.pass_to_HELAS(ind), coeff,
                                        self.write_obj(numerator.get_rep(ind))))
        #out.write('    mgDebug( 1, __FUNCTION__ );\n') # AV
        out.write('    return;\n') # AV
        ###return out.getvalue() # AV
        # AV check if one, two or half are used and need to be defined (ugly hack for #291: can this be done better?)
        out2 = StringIO()
        if 'one' in out.getvalue(): out2.write('    constexpr fptype one( 1. );\n')
        if 'two' in out.getvalue(): out2.write('    constexpr fptype two( 2. );\n')
        if 'half' in out.getvalue(): out2.write('    constexpr fptype half( 1. / 2. );\n')
        out2.write( out.getvalue() )
        return out2.getvalue()

    # AV - modify aloha_writers.WriteALOHA method (improve formatting)
    def write_MultVariable(self, obj, prefactor=True):
        """Turn a multvariable into a string"""
        mult_list = [self.write_variable_id(id) for id in obj]
        ###data = {'factors': '*'.join(mult_list)}
        data = {'factors': ' * '.join(mult_list)}
        if prefactor and obj.prefactor != 1:
            if obj.prefactor != -1:
                text = '%(prefactor)s * %(factors)s'
                data['prefactor'] = self.change_number_format(obj.prefactor)
            else:
                text = '-%(factors)s' # AV keep default (this is not used in eemumu)
        else:
            text = '%(factors)s'
        return text % data

    # AV - modify aloha_writers.WriteALOHA method (improve formatting)
    def write_MultContainer(self, obj, prefactor=True):
        """Turn a multvariable into a string"""
        mult_list = [self.write_obj(id) for id in obj]
        ###data = {'factors': '*'.join(mult_list)}
        data = {'factors': ' * '.join(mult_list)} # AV
        if prefactor and obj.prefactor != 1:
            if obj.prefactor != -1:
                text = '%(prefactor)s * %(factors)s'
                data['prefactor'] = self.change_number_format(obj.prefactor)
            else:
                text = '-%(factors)s' # AV keep default (this is not used in eemumu)
        else:
            text = '%(factors)s'
        return text % data

    # AV - new method (based on implementation of write_obj and write_MultVariable)
    def objIsSimpleVariable(self, obj) :
        ###print ( obj.vartype, obj.prefactor, len( obj ), obj ) # AV - FOR DEBUGGING
        return ( obj.vartype == 0 ) or ( obj.vartype == 2 and len( obj ) == 1 )

    # AV - modify aloha_writers.WriteALOHA method (improve formatting)
    # This affects 'V1[2] = ' and 'F1[2] = ' in HelAmps_sm.cc
    def write_obj_Add(self, obj, prefactor=True):
        """Turns addvariable into a string"""
        data = defaultdict(list)
        number = []
        [data[p.prefactor].append(p) if hasattr(p, 'prefactor') else number.append(p) for p in obj]
        file_str = StringIO()
        if prefactor and obj.prefactor != 1:
            formatted = self.change_number_format(obj.prefactor)
            if formatted.startswith(('+','-')):
                file_str.write('( %s )' % formatted)
            else:
                file_str.write(formatted)
            file_str.write(' * ( ')
        else:
            file_str.write('( ')
        ###print('."'+file_str.getvalue()+'"') # AV - FOR DEBUGGING
        first=True
        for value, obj_list in data.items():
            ###print('.."' + str(value) + '" "' + str(obj_list) + '"') # AV - FOR DEBUGGING
            add= ' + '
            if value not in  [-1,1]:
                nb_str = self.change_number_format(value)
                ###print('.>"' + nb_str + '"') # AV - FOR DEBUGGING
                if nb_str[0] in ['+', '-']:
                    if first: file_str.write(nb_str)
                    else : file_str.write(' ' + nb_str[0] + ' ' + nb_str[1:])
                elif first and nb_str == '( half * cI )':
                    file_str.write('half * cI')
                elif not first and nb_str == '( -half * cI )':
                    file_str.write(' - half * cI')
                else:
                    file_str.write('+' if first else ' + ')
                    file_str.write(nb_str)
                file_str.write(' * ')
                if len( obj_list ) > 1 or not self.objIsSimpleVariable( obj_list[0] ) : file_str.write('( ')
            elif value == -1:
                add = ' - '
                file_str.write('-' if first else ' - ')
            elif not first:
                file_str.write(' + ')
            else:
                file_str.write('')
            first = False
            # AV comment: write_obj here also adds calls declaration_add (via change_var_format) - example: OM3
            ###print('..."'+file_str.getvalue()+'"') # AV - FOR DEBUGGING
            file_str.write( add.join( [self.write_obj(obj, prefactor=False) for obj in obj_list] ) ) # NB: RECURSIVE! (write_obj_Add calls write_obj...)
            ###print('...."'+file_str.getvalue()+'"') # AV - FOR DEBUGGING
            if value not in [1,-1]:
                if len( obj_list ) > 1 or not self.objIsSimpleVariable( obj_list[0] ) : file_str.write(' )')
        if number:
            total = sum(number)
            file_str.write('+ %s' % self.change_number_format(total))
        file_str.write(' )')
        ###print('....."'+file_str.getvalue()+'"') # AV - FOR DEBUGGING
        return file_str.getvalue()

#------------------------------------------------------------------------------------

from os.path import join as pjoin

# AV - define a custom UFOModelConverter
# (NB: enable this via PLUGIN_ProcessExporter.create_model_class in output.py)
class PLUGIN_UFOModelConverter(PLUGIN_export_cpp.UFOModelConverterGPU):
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
    output_name = 'KOKKOS standalone'

    # AV - change defaults from export_cpp.UFOModelConverterGPU
    ###cc_ext = 'cu' # create HelAmps_sm.cu
    cc_ext = 'cc' # create HelAmps_sm.cc

    # AV - keep defaults from export_cpp.UFOModelConverterGPU
    ###cc_ext = 'cu'
    ###aloha_template_h = pjoin('gpu','cpp_hel_amps_h.inc')
    ###aloha_template_cc = pjoin('gpu','cpp_hel_amps_cc.inc')
    ###helas_h = pjoin('gpu', 'helas.h')
    helas_cc = pjoin('gpu', 'helas.cpp')

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
        res = res.replace('std::complex<','Kokkos::complex<') # custom simplex complex class (with constexpr arithmetics)
        if res == '' : res = '  // (none)'
        else : res = '  ' + res # add leading '  ' after the '// Model' line
        res = res.replace('\n','\n  ')
        res = res.replace(',',', ')
        return res

    # AV - overload export_cpp.UFOModelConverterCPP method (improve formatting)
    def write_set_parameters(self, params):
        res = super().write_set_parameters(params)
        res = res.replace('std::complex<','Kokkos::complex<') # custom simplex complex class (with constexpr arithmetics)
        if res == '' : res = '// (none)'
        res = res.replace('\n','\n  ')
        return res

    def write_hardcoded_parameters(self, params):
        pardef = super().write_parameters(params)
        parset = super().write_set_parameters(params)
        if ( pardef == '' ):
            assert( parset == '' ) # AV sanity check (both are empty)
            res = '// (none)\n'
            return res
        pardef = pardef.replace('std::complex<','Kokkos::complex<') # custom simplex complex class (with constexpr arithmetics)
        parset = parset.replace('std::complex<','Kokkos::complex<') # custom simplex complex class (with constexpr arithmetics)
        parset = parset.replace('sqrt(','constexpr_sqrt(') # constexpr sqrt (based on iterative Newton-Raphson approximation)
        parset = parset.replace('pow(','constexpr_pow(') # constexpr sqrt (based on iterative Newton-Raphson approximation)
        parset = parset.replace('(','( ')
        parset = parset.replace(')',' )')
        parset = parset.replace('+',' + ')
        parset = parset.replace('-',' - ')
        parset = parset.replace('e + ','e+') # fix exponents
        parset = parset.replace('e - ','e-') # fix exponents
        parset = parset.replace('=  + ','= +') # fix leading + in assignmments
        parset = parset.replace('=  - ','= -') # fix leading - in assignmments
        #parset = parset.replace('*',' * ')
        #parset = parset.replace('/',' / ')
        parset = parset.replace(',',', ')
        pardef_lines = {}
        for line in pardef.split('\n'):
            type, pars = line.rstrip(';').split(' ') # strip trailing ';'
            for par in pars.split(','):
                pardef_lines[par] = ( 'constexpr ' + type + ' ' + par )
        ###print( pardef_lines )
        parset_pars = []
        parset_lines = {}
        for line in parset.split('\n'):
            par, parval = line.split(' = ')
            if parval.startswith('slha.get_block_entry'): parval = parval.split(',')[2].lstrip(' ').rstrip(');') + ';'
            parset_pars.append( par )
            parset_lines[par] = parval # includes a trailing ';'
        ###print( parset_lines )
        assert( len(pardef_lines) == len(parset_lines) ) # AV sanity check (same number of parameters)
        res = '  '.join( pardef_lines[par] + ' = ' + parset_lines[par] + '\n' for par in parset_pars ) # no leading '  ' on first row
        res = res.replace(' ;',';')
        ###print(res); assert(False)
        return res

    def super_generate_parameters_class_files(self):
        """Create the content of the Parameters_model.h and .cc files"""
        replace_dict = self.default_replace_dict
        replace_dict['info_lines'] = PLUGIN_export_cpp.get_mg5_info_lines()
        replace_dict['model_name'] = self.model_name
        params_indep = [ line.replace('aS, ','')
                         for line in self.write_parameters(self.params_indep).split('\n') ]
        replace_dict['independent_parameters'] = '// Model parameters independent of aS\n  //double aS; // now retrieved event-by-event (as G) from Fortran (running alphas #373)\n' + '\n'.join( params_indep )
        replace_dict['independent_couplings'] = '// Model couplings independent of aS\n' + self.write_parameters(self.coups_indep)
        params_dep = [ '  //' + line[2:] + ' // now computed event-by-event (running alphas #373)' for line in self.write_parameters(self.params_dep).split('\n') ]
        replace_dict['dependent_parameters'] = '// Model parameters dependent on aS\n' + '\n'.join( params_dep )
        coups_dep = [ '  //' + line[2:] + ' // now computed event-by-event (running alphas #373)' for line in self.write_parameters(list(self.coups_dep.values())).split('\n') ]
        replace_dict['dependent_couplings'] = '// Model couplings dependent on aS\n' + '\n'.join( coups_dep )
        set_params_indep = [ line.replace('aS','//aS') + ' // now retrieved event-by-event (as G) from Fortran (running alphas #373)'
                             if line.startswith( '  aS =' ) else
                             line for line in self.write_set_parameters(self.params_indep).split('\n') ]
        replace_dict['set_independent_parameters'] = '\n'.join( set_params_indep )
        replace_dict['set_independent_couplings'] = self.write_set_parameters(self.coups_indep)
        replace_dict['set_dependent_parameters'] = self.write_set_parameters(self.params_dep)
        replace_dict['set_dependent_couplings'] = self.write_set_parameters(list(self.coups_dep.values()))
        print_params_indep = [ line.replace('std::cout','//std::cout') + ' // now retrieved event-by-event (as G) from Fortran (running alphas #373)'
                               if '"aS =' in line else
                               line for line in self.write_print_parameters(self.params_indep).split('\n') ]
        replace_dict['print_independent_parameters'] = '\n'.join( print_params_indep )
        replace_dict['print_independent_couplings'] = self.write_print_parameters(self.coups_indep)
        replace_dict['print_dependent_parameters'] = self.write_print_parameters(self.params_dep)
        replace_dict['print_dependent_couplings'] = self.write_print_parameters(list(self.coups_dep.values()))
        if 'include_prefix' not in replace_dict:
            replace_dict['include_prefix'] = ''
        hrd_params_indep = [ line.replace('constexpr','//constexpr') + ' // now retrieved event-by-event (as G) from Fortran (running alphas #373)' if 'aS =' in line else line for line in self.write_hardcoded_parameters(self.params_indep).split('\n') ]
        replace_dict['hardcoded_independent_parameters'] = '\n'.join( hrd_params_indep )
        replace_dict['hardcoded_independent_couplings'] = self.write_hardcoded_parameters(self.coups_indep)
        hrd_params_dep = [ line.replace('constexpr','//constexpr') + ' // now computed event-by-event (running alphas #373)' if line != '' else line for line in self.write_hardcoded_parameters(self.params_dep).split('\n') ]
        replace_dict['hardcoded_dependent_parameters'] = '\n'.join( hrd_params_dep )
        hrd_coups_dep = [ line.replace('constexpr','//constexpr') + ' // now computed event-by-event (running alphas #373)' if line != '' else line for line in self.write_hardcoded_parameters(list(self.coups_dep.values())).split('\n') ]
        replace_dict['hardcoded_dependent_couplings'] = '\n'.join( hrd_coups_dep )
        replace_dict['nicoup'] = len( self.coups_indep )
        if len( self.coups_indep ) > 0 :
            iicoup = [ '  constexpr size_t ixcoup_%s = %d + Parameters_%s_dependentCouplings::ndcoup; // out of ndcoup+nicoup' % (par.name, idx, self.model_name) for (idx, par) in enumerate(self.coups_indep) ]
            replace_dict['iicoup'] = '\n'.join( iicoup )
            icoupseticoup_hrdcod = [ '    ((FPType)Parameters_{2:s}::{0:s}.real(), (FPType)Parameters_{2:s}::{0:s}.imag()),'.format(par.name, idx, self.model_name) for (idx, par) in enumerate(self.coups_indep) ]
            replace_dict['icoupseticoup_hrdcod'] = '\n'.join( icoupseticoup_hrdcod )
        else:
            replace_dict['iicoup'] = '  // NB: there are no aS-independent couplings in this physics process'
            replace_dict['icoupseticoup_hrdcod'] = '    // NB: there are no aS-independent couplings in this physics process'
        replace_dict['ndcoup'] = len( self.coups_dep )
        if len( self.coups_dep ) > 0 :
            idcoup = [ '  constexpr size_t idcoup_%s = %d;' % (name, id) for (id, name) in enumerate(self.coups_dep) ]
            replace_dict['idcoup'] = '\n'.join( idcoup )
            dcoupdecl = [ '    cxtype %s;' % name for name in self.coups_dep ]
            replace_dict['dcoupdecl'] = '\n'.join( dcoupdecl )
            dcoupsetdpar = []
            foundG = False
            for line in self.write_hardcoded_parameters(self.params_dep).split('\n'):
                if line != '':
                    dcoupsetdpar.append( '    ' + line.replace('constexpr double', 'const FPType' if foundG else '//const FPType' ) )
                    if 'constexpr double G =' in line: foundG = True
            replace_dict['dcoupsetdpar'] = '  ' + '\n'.join( dcoupsetdpar )
            dcoupsetdcoup = [ '    ' + line.replace('constexpr Kokkos::complex<double> ','couplings[idcoup_')
                                           .replace(" = ","] = ")
                                           .replace('mdl_complexi', 'cI')
                                for line in self.write_hardcoded_parameters(list(self.coups_dep.values())).split('\n') if line != '' ]
            replace_dict['dcoupsetdcoup'] = '  ' + '\n'.join( dcoupsetdcoup )
        else:
            replace_dict['idcoup'] = '  // NB: there are no aS-dependent couplings in this physics process'
            replace_dict['dcoupdecl'] = '    // (none)'
            replace_dict['dcoupsetdpar'] = '      // (none)'
            replace_dict['dcoupsetdcoup'] = '      // (none)'
        if self.model_name == 'sm' :
            replace_dict['efterror'] = ''
        else:
            replace_dict['efterror'] = '\n#error This non-SM physics process only supports MGONGPU_HARDCODE_PARAM builds (#439): please run "make HRDCOD=1"'
        file_h = self.read_template_file(self.param_template_h) % replace_dict
        file_cc = self.read_template_file(self.param_template_cc) % replace_dict
        return file_h, file_cc

    # AV - overload export_cpp.UFOModelConverterCPP method (improve formatting)
    def generate_parameters_class_files(self):
        #file_h, file_cc = super().generate_parameters_class_files()
        file_h, file_cc = self.super_generate_parameters_class_files()
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
            res_strings.append('std::cout << std::setw( 20 ) << \"%s = \" << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << %s << std::endl;' % (param.name, param.name)) # AV
        if len(res_strings) == 0 : res_strings.append('// (none)')
        ##return '\n'.join(res_strings)
        return '\n  '.join(res_strings) # AV (why was this not necessary before?)

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
        replace_dict['info_lines'] = PLUGIN_export_cpp.get_mg5_info_lines()
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
        file_h_lines = file_h.split('\n')
        file_h = '\n'.join( file_h_lines) # skip the trailing '//---'
        file_h += file_cc # append the contents of HelAmps_sm.cc directly to HelAmps_sm.h!
        PLUGIN_writers.CPPWriter(model_h_file).writelines(file_h)
        logger.info('Created file %s in directory %s' \
                    % (os.path.split(model_h_file)[-1], os.path.split(model_h_file)[0] ) )

#------------------------------------------------------------------------------------

import madgraph.iolibs.files as files
import madgraph.various.misc as misc

# AV - define a custom OneProcessExporter
# (NB: enable this via PLUGIN_ProcessExporter.oneprocessclass in output.py)
# (NB: use this directly also in PLUGIN_UFOModelConverter.read_template_file)
# (NB: use this directly also in PLUGIN_GPUFOHelasCallWriter.super_get_matrix_element_calls)
class PLUGIN_OneProcessExporter(PLUGIN_export_cpp.OneProcessExporterGPU):
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

    # AV - overload export_cpp.OneProcessExporterGPU constructor (rename gCPPProcess to CPPProcess, set include_multi_channel)
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.process_class = 'CPPProcess'

    # AV - overload export_cpp.OneProcessExporterGPU method (indent comments in process_lines)
    def get_process_class_definitions(self, write=True):
        replace_dict = super().get_process_class_definitions(write=False)
        replace_dict['process_lines'] = replace_dict['process_lines'].replace('\n','\n  ')
        file = self.read_template_file(self.process_class_template) % replace_dict # HACK! ignore write=False case
        return file

    # AV - replace export_cpp.OneProcessExporterGPU method (fix gCPPProcess.cu)
    def get_process_function_definitions(self, write=True):
        """The complete class definition for the process"""
        replace_dict = super(PLUGIN_export_cpp.OneProcessExporterGPU,self).get_process_function_definitions(write=False)
        replace_dict['ncouplings'] = len(self.couplings2order)
        replace_dict['ncouplingstimes2'] = 2 * replace_dict['ncouplings']
        replace_dict['nparams'] = len(self.params2order)
        replace_dict['nmodels'] = replace_dict['nparams'] + replace_dict['ncouplings']
        replace_dict['coupling_list'] = ' '
        replace_dict['hel_amps_cc'] = '#include \"HelAmps_%s.cc\"' % self.model_name # AV
        coupling = [''] * len(self.couplings2order)
        params = [''] * len(self.params2order)
        for coup, pos in self.couplings2order.items():
            coupling[pos] = coup

        coupling_indep = [] # AV keep only the alphas-independent couplings #434
        for coup in coupling:
            keep = True
            # Use the same implementation as in UFOModelConverterCPP.prepare_couplings (assume self.model is the same)
            for key, coup_list in self.model['couplings'].items():
                if "aS" in key and coup in coup_list: keep = False
            if keep: coupling_indep.append( coup ) # AV only indep!
        
        ## NSN - Need access to independent couplings tIPC outside of CPPProcess for SYCL
        coup_str = ""
        #for i in range(len(self.couplings2order)):
        #    coup_str += "m_tIPC[%s] = cxmake( m_pars->%s );\n" % (i, coupling[i])
        if len(coupling_indep) > 0:
            for i in range(len(coupling_indep)):
                coup_str += "    m_tIPC[%s] = cxmake( m_pars->%s );\n" % (i, coupling_indep[i])
        else:
            coup_str = "    //m_tIPC[...] = ... ; // nicoup=0\n"

        self.number_dependent_couplings = len(coupling) - len(coupling_indep)
        self.number_independent_couplings = len(coupling_indep)

        for para, pos in self.params2order.items():
            params[pos] = para

        # NSN - Need access to tIPD outside of CPPProcess for SYCL
        param_str = ""
        if len(params) > 0:
            for i in range(len(self.params2order)):
                param_str += "    m_tIPD[%s] = (fptype)m_pars->%s;\n" % (i, params[i])
        else:
            parm_str += "    //m_tIPD[...] = ... ; // nparam=0\n"

        replace_dict['assign_coupling'] = coup_str + param_str
        file = self.read_template_file(self.process_definition_template) % replace_dict
        self.__process_function_definitions__ = file
        return file

    # AV - modify export_cpp.OneProcessExporterGPU method (add debug printouts for multichannel #342)
    def get_sigmaKin_lines(self, color_amplitudes, write=True):
        misc.sprint('Entering PLUGIN_OneProcessExporter.get_sigmaKin_lines')
        misc.sprint(self.include_multi_channel)
        misc.sprint(self.support_multichannel)
        return super().get_sigmaKin_lines(color_amplitudes, write)

    # AV - modify export_cpp.OneProcessExporterGPU method (fix gCPPProcess.cu)
    def get_all_sigmaKin_lines(self, color_amplitudes, class_name):
        """Get sigmaKin_process for all subprocesses for gCPPProcess.cu"""
        ret_lines = []
        if self.single_helicities:
            ret_lines.append("""
  // Evaluate |M|^2 for each subprocess
  // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event(s)
template <typename mom_t, typename ipc_t, typename ipd_t>
KOKKOS_INLINE_FUNCTION fptype calculate_wavefunctions(
  const mom_t& allmomenta,              // input: momenta
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
  fptype* allNumerators,                // output: multichannel numerators, running_sum_over_helicities
  fptype* allDenominators,              // output: multichannel denominators, running_sum_over_helicities
  const unsigned int channelId,         // input: multichannel channel id (1 to #diagrams); 0 to disable channel enhancement
#endif
  const short*  __restrict__ cHel,
  const ipc_t& COUPs,
  const ipd_t& cIPD
  )
{
  using namespace MG5_sm;
  fptype allMEs = 0;""")
            
            ret_lines.append("  // The number of colors")
            ret_lines.append("  constexpr int ncolor = %i;" % len(color_amplitudes[0]))
            ret_lines.append("""
  // Local TEMPORARY variables for a subset of Feynman diagrams in the given event (ievt)
  // [NB these variables are reused several times (and re-initialised each time) within the same event or event page]
  cxtype w[mgOnGpu::nwf][mgOnGpu::nw6]; // particle wavefunctions within Feynman diagrams (nw6 is often 6, the dimension of spin 1/2 or spin 1 particles)
  cxtype amp[1]; // invariant amplitude for one given Feynman diagram

  // Local variables for the given event (ievt)
  cxtype jamp[ncolor] = {}; // sum of the invariant amplitudes for all Feynman diagrams in the event or event page

  // === Calculate wavefunctions and amplitudes for all diagrams in all processes - Loop over nevt events ===
""")

            multi_channel = None
            if self.include_multi_channel:
                if not self.support_multichannel:
                    raise Exception("link with madevent not supported")
                multi_channel = self.get_multi_channel_dictionary(self.matrix_elements[0].get('diagrams'), self.include_multi_channel)
            helas_calls = self.helas_call_writer.get_matrix_element_calls(\
                                                    self.matrix_elements[0],
                                                    color_amplitudes[0],
                                                    multi_channel_map = multi_channel
                                                    )
            logger.debug("only one Matrix-element supported?")
            self.couplings2order = self.helas_call_writer.couplings2order
            self.params2order = self.helas_call_writer.params2order
            nwavefuncs = self.matrix_elements[0].get_number_of_wavefunctions()
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
                                                         class_name) for i, me in enumerate(self.matrix_elements)])
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
        """Generate mgOnGpuConfig.h, CPPProcess.cc, CPPProcess.h, check_sa.cc, gXXX.cu links"""
        misc.sprint('Entering PLUGIN_OneProcessExporter.generate_process_files')
        if self.include_multi_channel:
            misc.sprint('self.include_multi_channel is already defined: this is madevent+second_exporter mode')
        else:
            misc.sprint('self.include_multi_channel is not yet defined: this is standalone_kokkos mode') # see issue #473
        if self.matrix_elements[0].get('has_mirror_process'):
            self.matrix_elements[0].set('has_mirror_process', False)
            self.nprocesses/=2

        ## NSN Explicity copy method and generate cc file first (need params2order for h file)
        #super(export_cpp.OneProcessExporterGPU, self).generate_process_files()

        #"""Generate the .h and .cc files needed for C++, for the
        #processes described by multi_matrix_element"""

        # Create the files
        if not os.path.isdir(os.path.join(self.path, self.process_dir)):
            os.makedirs(os.path.join(self.path, self.process_dir))
        filename = os.path.join(self.path, self.process_dir,
                                '%s.%s' % (self.process_class, self.cc_ext))

        self.write_process_cc_file(PLUGIN_writers.CPPWriter(filename))


        if not os.path.isdir(os.path.join(self.path, self.include_dir)):
            os.makedirs(os.path.join(self.path, self.include_dir))
        filename = os.path.join(self.path, self.include_dir,
                                '%s.h' % self.process_class)

        self.write_process_h_file(PLUGIN_writers.CPPWriter(filename))


        logger.info('Created files %(process)s.h and %(process)s.cc in' % \
                    {'process': self.process_class} + \
                    ' directory %(dir)s' % {'dir': os.path.split(filename)[0]})

        self.edit_CMakeLists()
        self.edit_check_sa()
        self.edit_mgonGPU()
        #self.edit_processidfile() # AV new file (NB this is Sigma-specific, should not be a symlink to Subprocesses)
        #self.edit_testxxx() # AV new file (NB this is generic in Subprocesses and then linked in Sigma-specific)
        # Add symbolic links
        # NB: symlink of kokkos.mk to makefile is overwritten by madevent makefile if this exists (#480)
        # NB: this relies on the assumption that kokkos code is generated before madevent code
        files.ln(pjoin(self.path, 'kokkos.mk'), self.path, 'makefile')
        #files.ln(pjoin(self.path, 'check_sa.cc'), self.path, 'gcheck_sa.cu')
        #files.ln(pjoin(self.path, 'CPPProcess.cc'), self.path, 'gCPPProcess.cu')
        #files.ln(pjoin(self.path, 'CrossSectionKernels.cc'), self.path, 'gCrossSectionKernels.cu')
        #files.ln(pjoin(self.path, 'MatrixElementKernels.cc'), self.path, 'gMatrixElementKernels.cu')
        #files.ln(pjoin(self.path, 'RamboSamplingKernels.cc'), self.path, 'gRamboSamplingKernels.cu')
        #files.ln(pjoin(self.path, 'RandomNumberKernels.cc'), self.path, 'gRandomNumberKernels.cu')
        #files.ln(pjoin(self.path, 'BridgeKernels.cc'), self.path, 'gBridgeKernels.cu')

    # SR - generate CMakeLists.txt file inside the P* directory
    def edit_CMakeLists(self):
        """Generate CMakeLists.txt"""
        misc.sprint('Entering PLUGIN_OneProcessExporter.edit_CMakeLists')
        template = open(pjoin(self.template_path,'CMake/SubProcesses/CMakeLists_P.txt'),'r').read()
        ff = open(pjoin(self.path, 'CMakeLists.txt'),'w')
        ff.write(template)
        ff.close()

    # AV - replace the export_cpp.OneProcessExporterGPU method (invert .cc/.cu, add debug printouts)
    def edit_check_sa(self):
        """Generate check_sa.cc and fcheck_sa.f"""
        misc.sprint('Entering PLUGIN_OneProcessExporter.edit_check_sa')
        ff = open(pjoin(self.path, 'check_sa.cc'),'w')
        template = open(pjoin(self.template_path,'gpu','check.cpp'),'r').read()
        ff.write(template)
        ff.close()

        replace_dict = {}
        replace_dict['nexternal'], _ = self.matrix_elements[0].get_nexternal_ninitial()
        # replace_dict['model'] = self.model_name
        # replace_dict['numproc'] = len(self.matrix_elements)
        ff = open(pjoin(self.path, 'fcheck_sa.f'),'w')
        template = open(pjoin(self.template_path,'gpu','fcheck_sa.f'),'r').read()
        ff.write(template % replace_dict)
        ff.close()

    # AV - replace the export_cpp.OneProcessExporterGPU method (add debug printouts and multichannel handling #473) 
    def edit_mgonGPU(self):
        """Generate mgOnGpuConfig.h"""
        misc.sprint('Entering PLUGIN_OneProcessExporter.edit_mgonGPU')
        template = open(pjoin(self.template_path,'gpu','mgOnGpuConfig.h'),'r').read()
        replace_dict = {}
        nexternal, nincoming = self.matrix_elements[0].get_nexternal_ninitial()
        replace_dict['number_dependent_couplings'] = self.number_dependent_couplings
        replace_dict['number_independent_couplings'] = self.number_independent_couplings
        replace_dict['nincoming'] = nincoming
        replace_dict['noutcoming'] = nexternal - nincoming

        # Number of helicity combinations
        replace_dict['nbhel'] = \
                            self.matrix_elements[0].get_helicity_combinations()
        replace_dict['nwavefunc'] = \
                          self.matrix_elements[0].get_number_of_wavefunctions()
        replace_dict['wavefuncsize'] = 6

        replace_dict['ncouplings'] = len(self.couplings2order)
        replace_dict['ncouplingstimes2'] = 2 * replace_dict['ncouplings']
        replace_dict['nparams'] = len(self.params2order)
        
        
        if self.include_multi_channel:
            replace_dict['mgongpu_supports_multichannel'] = '#define MGONGPU_SUPPORTS_MULTICHANNEL 1'
        else:
            replace_dict['mgongpu_supports_multichannel'] = '#undef MGONGPU_SUPPORTS_MULTICHANNEL'

        ff = open(pjoin(self.path, '..','..','src','mgOnGpuConfig.h'),'w')
        ff.write(template % replace_dict)
        ff.close()

    # AV - new method
    def edit_processidfile(self):
        """Generate epoch_process_id.h"""
        misc.sprint('Entering PLUGIN_OneProcessExporter.edit_processidfile')
        template = open(pjoin(self.template_path,'gpu','epoch_process_id.h'),'r').read()
        replace_dict = {}
        replace_dict['processid'] = self.get_process_name()
        replace_dict['processid_uppercase'] = self.get_process_name().upper()
        ff = open(pjoin(self.path, 'epoch_process_id.h'),'w')
        ff.write(template % replace_dict)
        ff.close()

    # AV - new method
    def edit_testxxx(self):
        """Generate testxxx.cc"""
        misc.sprint('Entering PLUGIN_OneProcessExporter.edit_testxxx')
        template = open(pjoin(self.template_path,'gpu','testxxx.cc'),'r').read()
        replace_dict = {}
        replace_dict['model_name'] = self.model_name
        ff = open(pjoin(self.path, '..', 'testxxx.cc'),'w')
        ff.write(template % replace_dict)
        ff.close()

    # AV - overload the export_cpp.OneProcessExporterGPU method (add debug printout and truncate last \n)
    # [*NB export_cpp.UFOModelConverterGPU.write_process_h_file is not called!*]
    def write_process_h_file(self, writer):
        """Generate final CPPProcess.h"""
        misc.sprint('Entering PLUGIN_OneProcessExporter.write_process_h_file')
        replace_dict = super(PLUGIN_export_cpp.OneProcessExporterGPU, self).write_process_h_file(False)
        #replace_dict2 = super(PLUGIN_export_cpp.OneProcessExporterGPU,self).get_process_function_definitions(write=False)
        replace_dict['helamps_h'] = "\n#include \"HelAmps_%s.h\"" % self.model_name

        # Kokkos puts source code in header and the helicities are set in the process_function_definitions
        #cc_replace_dict = super(PLUGIN_export_cpp.OneProcessExporterGPU, self).write_process_cc_file(False)
        #replace_dict['process_function_definitions'] = cc_replace_dict['process_function_definitions'] 
        replace_dict['process_function_definitions'] = self.__process_function_definitions__

        #Set helicities
        replace_dict['all_helicities'] = self.get_helicity_matrix(self.matrix_elements[0])
        # replace_dict['all_helicities'] = replace_dict['all_helicities'] .replace("{", "")
        # replace_dict['all_helicities'] = replace_dict['all_helicities'] .replace("}", "")
        # replace_dict['all_helicities'] = replace_dict['all_helicities'] .replace("helicities", "constexpr T helicities[mgOnGpu::ncomb][mgOnGpu::npar] {")
        # replace_dict['all_helicities'] = replace_dict['all_helicities'] .replace("    ", "  ")
        # replace_dict['all_helicities'] = replace_dict['all_helicities'] .replace(", constexpr", "constexpr")
        # replace_dict['all_helicities'] = replace_dict['all_helicities'] .replace(";", "")

        #Set hardcoded parameters
        params = [''] * len(self.params2order)
        for para, pos in self.params2order.items():
            params[pos] = para

        misc.sprint(replace_dict)
        if writer:
            file = self.read_template_file(self.process_template_h) % replace_dict
            # Write the file
            writer.writelines(file)
        else:
            return replace_dict
        #out = super().write_process_h_file(writer)
        writer.seek(-1, os.SEEK_CUR)
        writer.truncate()

    # AV - replace the export_cpp.OneProcessExporterGPU method (replace HelAmps.cu by HelAmps.cc)
    def super_write_process_cc_file(self, writer):
        """Write the class member definition (.cc) file for the process described by matrix_element"""
        replace_dict = super(PLUGIN_export_cpp.OneProcessExporterGPU, self).write_process_cc_file(False)
        ###replace_dict['hel_amps_def'] = '\n#include \"../../src/HelAmps_%s.cu\"' % self.model_name
        replace_dict['hel_amps_h'] = '#include \"HelAmps_%s.h\"' % self.model_name # AV
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

    # AV - replace the export_cpp.OneProcessExporterGPU method (improve formatting? actually keep all defaults!)
    # [NB this is used in uu~>tt~ but not in gg>tt~ or e+e->mu+mu-, see issue #337]
    @staticmethod
    def coeff(ff_number, frac, is_imaginary, Nc_power, Nc_value=3):
        """Returns a nicely formatted string for the coefficients in JAMP lines"""
        total_coeff = ff_number * frac * Fraction(Nc_value) ** Nc_power
        if total_coeff == 1:
            if is_imaginary:
                return '+cxtype(0,1)*' # AV keep default (this is not used in eemumu - should use cI eventually)
            else:
                return '+' # AV keep default (this is not used in eemumu)
        elif total_coeff == -1:
            if is_imaginary:
                return '-cxtype(0,1)*' # AV keep default (this is not used in eemumu - should use cI eventually)
            else:
                return '-' # AV keep default (eg jamp[0] += -amp[0])
        ###assert(False) # [this had been inserted to check if coeff is used at all, it is used in uu~>tt~, see #337]
        res_str = '%+i.' % total_coeff.numerator
        if total_coeff.denominator != 1:
            # Check if total_coeff is an integer
            res_str = res_str + '/%i.' % total_coeff.denominator
        if is_imaginary:
            res_str = res_str + '*cxtype(0,1)'
        return res_str + '*' # AV keep default (this is not used in eemumu)

    # AV - replace the export_cpp.OneProcessExporterCPP method (fix fptype and improve formatting)
    def get_color_matrix_lines(self, matrix_element):
        """Return the color matrix definition lines for this matrix element. Split rows in chunks of size n."""
        import madgraph.core.color_algebra as color
        if not matrix_element.get('color_matrix'):
            return '\n'.join(['      static constexpr fptype denom[1] = {1.};', 'static const fptype cf[1][1] = {1.};'])
        else:
            color_denominators = matrix_element.get('color_matrix').\
                                                 get_line_denominators()
            denom_string = '      static constexpr fptype denom[ncolor] = { %s }; // 1-D array[%i]' \
                           % ( ', '.join(['%i' % denom for denom in color_denominators]), len(color_denominators) )
            matrix_strings = []
            my_cs = color.ColorString()
            for index, denominator in enumerate(color_denominators):
                # Then write the numerators for the matrix elements
                num_list = matrix_element.get('color_matrix').get_line_numerators(index, denominator)
                matrix_strings.append('{ %s }' % ', '.join(['%d' % i for i in num_list]))
            matrix_string = '      static constexpr fptype cf[ncolor][ncolor] = '
            if len( matrix_strings ) > 1 : matrix_string += '{\n        ' + ',\n        '.join(matrix_strings) + ' };'
            else: matrix_string += '{ ' + matrix_strings[0] + ' };'
            matrix_string += ' // 2-D array[%i][%i]' % ( len(color_denominators), len(color_denominators) )
            denom_comment = '\n      // The color denominators (initialize all array elements, with ncolor=%i)\n      // [NB do keep \'static\' for these constexpr arrays, see issue #283]\n' % len(color_denominators)
            matrix_comment = '\n      // The color matrix (initialize all array elements, with ncolor=%i)\n      // [NB do keep \'static\' for these constexpr arrays, see issue #283]\n' % len(color_denominators)
            denom_string = denom_comment + denom_string
            matrix_string = matrix_comment + matrix_string
            return '\n'.join([denom_string, matrix_string])

    # AV - replace the export_cpp.OneProcessExporterGPU method (improve formatting)
    def get_initProc_lines(self, matrix_element, color_amplitudes):
        """Get initProc_lines for function definition for gCPPProcess::initProc"""
        initProc_lines = []
        initProc_lines.append("// Set external particle masses for this matrix element")
        for part in matrix_element.get_external_wavefunctions():
            ###initProc_lines.append("mME.push_back(pars->%s);" % part.get('mass'))
            initProc_lines.append("    m_masses.push_back( m_pars->%s );" % part.get('mass')) # AV
        # initProc_lines.append('  Kokkos::deep_copy(cmME,hmME);')
        ###for i, colamp in enumerate(color_amplitudes):
        ###    initProc_lines.append("jamp2[%d] = new double[%d];" % (i, len(colamp))) # AV - this was commented out already
        return "\n".join(initProc_lines)

    # AV - replace the export_cpp.OneProcessExporterCPP method (improve formatting)
    def get_helicity_matrix(self, matrix_element):
        """Return the Helicity matrix definition lines for this matrix element"""
        ###helicity_line = "static const int helicities[ncomb][nexternal] = {";
        ###helicity_line = "    constexpr short helicities[ncomb][mgOnGpu::npar] = {\n      "; # AV (this is tHel)
        #helicity_line = "template <typename T>\nconstexpr T helicities[mgOnGpu::ncomb][mgOnGpu::npar] {\n  "; # NSN SYCL needs access to tHel outside CPPProcess
        helicity_line = ""; # NSN SYCL needs access to tHel outside CPPProcess
        helicity_line_list = []
        for helicities in matrix_element.get_helicity_matrix(allow_reverse=False):
            #helicity_line_list.append( "{" + ", ".join(['%d'] * len(helicities)) % tuple(helicities) + "}" ) # AV"
            helicity_line_list.append( "    " + ", ".join(['%d'] * len(helicities)) % tuple(helicities) )
        #return helicity_line + ",\n  ".join(helicity_line_list) + "\n};" # AV
        return helicity_line + ",\n  ".join(helicity_line_list)

    # AV - overload the export_cpp.OneProcessExporterGPU method (just to add some comments...)
    def get_reset_jamp_lines(self, color_amplitudes):
        """Get lines to reset jamps"""
        ret_lines = super().get_reset_jamp_lines(color_amplitudes)
        if ret_lines != '' : ret_lines = '    // Reset jamp (reset color flows)\n' + ret_lines # AV THIS SHOULD NEVER HAPPEN!
        return ret_lines

#------------------------------------------------------------------------------------

import madgraph.core.helas_objects as helas_objects
import madgraph.iolibs.helas_call_writers as helas_call_writers

# AV - define a custom HelasCallWriter
# (NB: enable this via PLUGIN_ProcessExporter.helas_exporter in output.py - this fixes #341)
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
    def format_coupling(self, call):
        """Format the coupling so any minus signs are put in front"""
        import re
        ###print(call) # FOR DEBUGGING
        model = self.get('model')
        if not hasattr(self, 'couplings2order'):
            self.couplings2order = {}
            self.params2order = {}
        for coup in re.findall(self.findcoupling, call):
            if coup == 'ZERO':
                ###call = call.replace('pars->ZERO', '0.')
                call = call.replace('m_pars->ZERO', '0.') # AV
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
                name = 'cIPD'
            else:
                alias = self.couplings2order
                name = 'cIPC'
            if coup not in alias:
                alias[coup] = len(alias)
            if name == 'cIPD':
                call = call.replace('m_pars->%s%s' % (sign, coup),
                                    '%s%s[%s]' % (sign, name, alias[coup]))
            else:
                call = call.replace('m_pars->%s%s' % (sign, coup),
                                    '%sCOUPs[%s]' % (sign, alias[coup]))
        return call

    # AV - new method for formatting wavefunction/amplitude calls
    # [It would be too complex to modify them in helas_objects.HelasWavefunction/Amplitude.get_call_key]
    @staticmethod
    def format_call(call):
        return call.replace('(','( ').replace(')',' )').replace(',',', ')

    # AV - replace helas_call_writers.GPUFOHelasCallWriter method (improve formatting)
    def super_get_matrix_element_calls(self, matrix_element, color_amplitudes, multi_channel_map=False):
        """Return a list of strings, corresponding to the Helas calls for the matrix element"""
        import madgraph.core.helas_objects as helas_objects
        import madgraph.loop.loop_helas_objects as loop_helas_objects
        assert isinstance(matrix_element, helas_objects.HelasMatrixElement), \
               '%s not valid argument for get_matrix_element_calls' % \
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
        res.append('// Reset color flows (reset jamp) at the beginning of a new event or event page')
        res.append('for( int i=0; i<ncolor; i++ ){ jamp[i] = cxmake( 0, 0 ); }')
        diagrams = matrix_element.get('diagrams')
        diag_to_config = {}
        if multi_channel_map:
            for config in sorted(multi_channel_map.keys()):
                amp = [a.get('number') for a in \
                                  sum([diagrams[idiag].get('amplitudes') for \
                                       idiag in multi_channel_map[config]], [])]
                diag_to_config[amp[0]] = config
        id_amp = 0
        for diagram in matrix_element.get('diagrams'):
            ###print('DIAGRAM %3d: #wavefunctions=%3d, #diagrams=%3d' %
            ###      (diagram.get('number'), len(diagram.get('wavefunctions')), len(diagram.get('amplitudes')) )) # AV - FOR DEBUGGING
            res.append('\n  // *** DIAGRAM %d OF %d ***' % (diagram.get('number'), len(matrix_element.get('diagrams'))) ) # AV
            res.append('\n  // Wavefunction(s) for diagram number %d' % diagram.get('number')) # AV
            res.extend([ self.get_wavefunction_call(wf) for wf in diagram.get('wavefunctions') ]) # AV new: avoid format_call
            if len(diagram.get('wavefunctions')) == 0 : res.append('// (none)') # AV
            if res[-1][-1] == '\n' : res[-1] = res[-1][:-1]
            res.append('\n  // Amplitude(s) for diagram number %d' % diagram.get('number'))
            for amplitude in diagram.get('amplitudes'):
                id_amp +=1
                namp = amplitude.get('number')
                amplitude.set('number', 1)
                res.append(self.get_amplitude_call(amplitude)) # AV new: avoid format_call
                if multi_channel_map: # different code bases #473 (assume this is the same as self.include_multi_channel...)
                    if id_amp in diag_to_config:
                        res.append("#ifdef MGONGPU_SUPPORTS_MULTICHANNEL")
                        res.append("if( channelId == %i ) allNumerators[0] += cxabs2( amp[0] );" % diagram.get('number'))
                        res.append("if( channelId != 0 ) allDenominators[0] += cxabs2( amp[0] );")
                        res.append("#endif")
                else:
                    res.append("#ifdef MGONGPU_SUPPORTS_MULTICHANNEL")
                    res.append("// Here the code base generated with multichannel support updates numerators and denominators (#473)")
                    res.append("#endif")
                for njamp, coeff in color[namp].items():
                    scoeff = PLUGIN_OneProcessExporter.coeff(*coeff) # AV
                    if scoeff[0] == '+' : scoeff = scoeff[1:]
                    scoeff = scoeff.replace('(','( ')
                    scoeff = scoeff.replace(')',' )')
                    scoeff = scoeff.replace(',',', ')
                    scoeff = scoeff.replace('*',' * ')
                    scoeff = scoeff.replace('/',' / ')
                    if scoeff.startswith('-'): res.append("jamp[%s] -= %samp[0];" % (njamp, scoeff[1:])) # AV
                    else: res.append("jamp[%s] += %samp[0];" % (njamp, scoeff)) # AV
            if len(diagram.get('amplitudes')) == 0 : res.append('// (none)') # AV
        ###res.append('\n    // *** END OF DIAGRAMS ***' ) # AV - no longer needed ('COLOR ALGEBRA BELOW')
        return res

    # AV - overload helas_call_writers.GPUFOHelasCallWriter method (improve formatting)
    def get_matrix_element_calls(self, matrix_element, color_amplitudes, multi_channel_map=False):
        """Return a list of strings, corresponding to the Helas calls for the matrix element"""
        res = self.super_get_matrix_element_calls(matrix_element, color_amplitudes, multi_channel_map)
        for i, item in enumerate(res):
            ###print(item) # FOR DEBUGGING
            if item.startswith('# Amplitude'): item='//'+item[1:] # AV replace '# Amplitude' by '// Amplitude'
            if not item.startswith('\n') and not item.startswith('#'): res[i]='  '+item
        return res

    # AV - replace helas_call_writers.GPUFOHelasCallWriter method (improve formatting)
    # [GPUFOHelasCallWriter.format_coupling is called by GPUFOHelasCallWriter.get_external_line/generate_helas_call]
    # [GPUFOHelasCallWriter.get_external_line is called by GPUFOHelasCallWriter.get_external]
    # [=> GPUFOHelasCallWriter.get_external is called by GPUFOHelasCallWriter.generate_helas_call]
    # [GPUFOHelasCallWriter.generate_helas_call is called by UFOHelasCallWriter.get_wavefunction_call/get_amplitude_call]
    first_get_external = True
    def get_external(self, wf, argument):
        line = self.get_external_line(wf, argument)
        split_line = line.split(',')
        split_line = [ str.lstrip(' ').rstrip(' ') for str in split_line] # AV
        # (AV join using ',': no need to add a space as this is done by format_call later on)
        line = ', '.join(split_line)
        if self.first_get_external and ( ( 'mzxxx' in line ) or ( 'pzxxx' in line ) or ( 'xzxxx' in line ) ) :
            self.first_get_external = False
            line2 = line.replace('mzxxx','xxxxx').replace('pzxxx','xxxxx').replace('xzxxx','xxxxx')
            line2 = line2[:line2.find('// NB')]
            split_line2 = line2.split(',')
            split_line2 = [ str.lstrip(' ').rstrip(' ') for str in split_line2] # AV
            split_line2.insert(1, '0') # add parameter fmass=0
            line2 = ', '.join(split_line2)
            text = '#if not defined KOKKOS_ENABLE_CUDA\n      %s\n#else\n      if ( ievt %% 2 == 0 )\n        %s\n      else\n        %s\n#endif\n'
            return text % (line, line, line2)
        text = '%s\n'
        return text % line

    # AV - replace helas_call_writers.GPUFOHelasCallWriter method (vectorize w)
    # This is the method that creates the ixxx/oxxx function calls in calculate_wavefunctions
    # [GPUFOHelasCallWriter.get_external_line is called by GPUFOHelasCallWriter.get_external]
    # [GPUFOHelasCallWriter.get_external (adding #ifdef CUDA) is called by GPUFOHelasCallWriter.generate_helas_call]
    # [GPUFOHelasCallWriter.generate_helas_call is called by UFOHelasCallWriter.get_wavefunction_call/get_amplitude_call]
    def get_external_line(self, wf, argument):
        call = ''
        call = call + helas_call_writers.HelasCallWriter.mother_dict[\
                argument.get_spin_state_number()].lower() 
        if wf.get('mass').lower() != 'zero' or argument.get('spin') != 2:
            # Fill out with X up to 6 positions
            call = call + 'x' * (6 - len(call))
            # Specify namespace for Helas calls
            call = call + "(Kokkos::subview(allmomenta,%d,Kokkos::ALL),"
            if argument.get('spin') != 1:
                # For non-scalars, need mass and helicity
                call = call + "m_pars->%s, cHel[%d],"
            else:
                # AV This seems to be for scalars (spin==1???), pass neither mass nor helicity (#351)
                ###call = call + 'm_pars->%s,'
                call = call
            call = call + '%+d, w[%d]);' # AV vectorize
            if argument.get('spin') == 1:
                # AV This seems to be for scalars (spin==1???), pass neither mass nor helicity (#351)
                return call % \
                                (
                                 wf.get('number_external')-1,
                                 wf.get('mass'),
                                 wf.get('number_external')-1,
                                 # For boson, need initial/final here
                                 (-1) ** (wf.get('state') == 'initial'),
                                 wf.get('me_id')-1,
                                 wf.get('number_external')-1)
            elif argument.is_boson():
                misc.sprint(call)
                misc.sprint( (
                                 wf.get('number_external')-1,
                                 wf.get('mass'),
                                 wf.get('number_external')-1,
                                 # For boson, need initial/final here
                                 (-1) ** (wf.get('state') == 'initial'),
                                 wf.get('me_id')-1))
                return  self.format_coupling(call % \
                                (wf.get('number_external')-1,
                                 wf.get('mass'),
                                 wf.get('number_external')-1,
                                 # For boson, need initial/final here
                                 (-1) ** (wf.get('state') == 'initial'),
                                 wf.get('me_id')-1))
            else:
                return self.format_coupling(call % \
                                (wf.get('number_external')-1,
                                 wf.get('mass'),
                                 wf.get('number_external')-1,
                                 # For fermions, need particle/antiparticle
                                 - (-1) ** wf.get_with_flow('is_part'),
                                 wf.get('me_id')-1))
        else:
            if wf.get('number_external') == 1:
                call += 'pz'
            elif wf.get('number_external') == 2:
                call += 'mz'
            else:
                call += 'xz'
                comment = '' # AV
            call = call + 'x' * (6 - len(call))
            if wf.get('number_external') == 1 or wf.get('number_external') == 2: # AV
                comment = ' // NB: ' + call + ' only uses pz' # AV skip '(not E,px,py)' to avoid interference with comma parsing in get_external
            # Specify namespace for Helas calls
            call = call + '( Kokkos::subview(allmomenta,%d,Kokkos::ALL), cHel[%d], %+d, w[%d] );' + comment # AV vectorize and add comment
            return self.format_coupling(call % \
                                (wf.get('number_external')-1,
                                 wf.get('number_external')-1,
                                 # For fermions, need particle/antiparticle
                                 - (-1) ** wf.get_with_flow('is_part'),
                                 wf.get('me_id')-1))

    # AV - replace helas_call_writers.GPUFOHelasCallWriter method (vectorize w and amp)
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
        if not isinstance(argument, helas_objects.HelasWavefunction) and \
           not isinstance(argument, helas_objects.HelasAmplitude):
            raise self.PhysicsObjectError('get_helas_call must be called with wavefunction or amplitude')
        call = ''
        call_function = None
        if isinstance(argument, helas_objects.HelasAmplitude) and \
           argument.get('interaction_id') == 0:
            call = '#'
            call_function = lambda amp: call
            self.add_amplitude(argument.get_call_key(), call_function)
            return
        if isinstance(argument, helas_objects.HelasWavefunction) and \
               not argument.get('mothers'):
            # String is just ixxxxx, oxxxxx, vxxxxx or sxxxxx
            call_function = lambda wf: self.get_external(wf, argument)
        else:
            if isinstance(argument, helas_objects.HelasWavefunction):
                outgoing = argument.find_outgoing_number()
            else:
                outgoing = 0
            # Check if we need to append a charge conjugation flag
            l = [str(l) for l in argument.get('lorentz')]
            flag = [] 
            if argument.needs_hermitian_conjugate():
                flag = ['C%d' % i for i in argument.get_conjugate_index()]
            # Creating line formatting:
            # (AV NB: in the default code these two branches were identical, use a single branch)
            ###if isinstance(argument, helas_objects.HelasWavefunction): # AV e.g. FFV1P0_3 (output is wavefunction)
            ###    call = '%(routine_name)s(%(wf)s%(coup)s%(mass)s%(out)s);'
            ###else: # AV e.g. FFV1_0 (output is amplitude)
            ###    call = '%(routine_name)s(%(wf)s%(coup)s%(mass)s%(out)s);'
            call = '%(routine_name)s( %(wf)s%(coup)s%(mass)s%(out)s );'
            # compute wf
            arg = {'routine_name': aloha_writers.combine_name('%s' % l[0], l[1:], outgoing, flag, True),
                   'wf': ("w[%%(%d)d], " * len(argument.get('mothers'))) % tuple(range(len(argument.get('mothers')))), # AV
                   'coup': ("m_pars->%%(coup%d)s, " * len(argument.get('coupling'))) % tuple(range(len(argument.get('coupling')))) # AV
                   }
            if isinstance(argument, helas_objects.HelasWavefunction):
                ###arg['out'] = 'w[%(out)d]'
                arg['out'] = 'w[%(out)d]'
                if aloha.complex_mass:
                    arg['mass'] = "m_pars->%(CM)s, "
                else:
                    arg['mass'] = "m_pars->%(M)s, m_pars->%(W)s, "
            else:        
                arg['out'] = '&amp[%(out)d]'
                arg['out2'] = 'amp[%(out)d]'
                arg['mass'] = ''
            call = call % arg
            # Now we have a line correctly formatted
            call_function = lambda wf: self.format_coupling(
                                         call % wf.get_helas_call_dict(index=0))
        # Add the constructed function to wavefunction or amplitude dictionary
        if isinstance(argument, helas_objects.HelasWavefunction):
            self.add_wavefunction(argument.get_call_key(), call_function)
        else:
            self.add_amplitude(argument.get_call_key(), call_function)

# AV - use the custom HelasCallWriter
DEFAULT_GPUFOHelasCallWriter = helas_call_writers.GPUFOHelasCallWriter
helas_call_writers.GPUFOHelasCallWriter = PLUGIN_GPUFOHelasCallWriter

#------------------------------------------------------------------------------------
