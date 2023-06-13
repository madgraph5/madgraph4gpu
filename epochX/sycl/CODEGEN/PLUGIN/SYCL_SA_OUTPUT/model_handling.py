# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: O. Mattelaer (Sep 2021) for the MG5aMC CUDACPP plugin.
# Further modified by: O. Mattelaer, A. Valassi (2021-2023) for the MG5aMC CUDACPP plugin.
#
# Copyright (C) 2021-2023 Argonne National Laboratory.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Modified by: N. Nichols (2021-2023) for the MG5aMC SYCL plugin.

import os
pjoin = os.path.join

# Use templates for source code, scripts and Makefiles from PLUGINDIR instead of MG5DIR
PLUGINDIR = os.path.dirname( __file__ )

# Create a plugin-specific logger
import logging
logger = logging.getLogger('madgraph.PLUGIN.SYCL_SA_OUTPUT.model_handling')

#------------------------------------------------------------------------------------

# Import the independent 2nd copy of the export_cpp module (as PLUGIN_export_cpp), previously loaded in output.py
#import madgraph.iolibs.export_cpp as export_cpp # 1st copy
##import madgraph.iolibs.export_cpp as PLUGIN_export_cpp # this is not enough to define an independent 2nd copy: id(export_cpp)==id(PLUGIN_export_cpp)
import PLUGIN.SYCL_SA_OUTPUT.PLUGIN_export_cpp as PLUGIN_export_cpp # 2nd copy loaded in the plugin's output.py
#print('id(export_cpp)=%s'%id(export_cpp))
#print('id(PLUGIN_export_cpp)=%s'%id(PLUGIN_export_cpp))

#------------------------------------------------------------------------------------

# Modify export_cpp.get_mg5_info_lines (replace '# ' by '//')

def PLUGIN_get_mg5_info_lines():
    return DEFAULT_get_mg5_info_lines().replace('# ','//')

DEFAULT_get_mg5_info_lines = PLUGIN_export_cpp.get_mg5_info_lines
PLUGIN_export_cpp.get_mg5_info_lines = PLUGIN_get_mg5_info_lines

#------------------------------------------------------------------------------------

# Load an independent 2nd copy of the writers module (as PLUGIN_writers) and use that within the plugin (workaround for #341)
# See https://stackoverflow.com/a/11285504
#import madgraph.iolibs.file_writers as writers # 1st copy
import sys
import importlib.util
SPEC_WRITERS = importlib.util.find_spec('madgraph.iolibs.file_writers')
PLUGIN_writers = importlib.util.module_from_spec(SPEC_WRITERS)
SPEC_WRITERS.loader.exec_module(PLUGIN_writers)
#sys.modules['PLUGIN.SYCL_SA_OUTPUT.PLUGIN_writers'] = PLUGIN_writers # would allow 'import PLUGIN.SYCL_SA_OUTPUT.PLUGIN_writers' (not needed)
del SPEC_WRITERS

# Use the independent 2nd copy of the writers module within the PLUGIN_export_cpp module (workaround for #341)
###DEFAULT_writers = PLUGIN_export_cpp.writers # not needed
PLUGIN_export_cpp.writers = PLUGIN_writers

#------------------------------------------------------------------------------------

# Modify writers.FileWriter.__init__ (add a debug printout)
def PLUGIN_FileWriter__init__( self, name, opt = 'w' ):
    print( 'FileWriter %s for %s'%( type(self), name) )
    return DEFAULT_FileWriter__init__( self, name, opt )

DEFAULT_FileWriter__init__ = PLUGIN_writers.FileWriter.__init__
PLUGIN_writers.FileWriter.__init__ = PLUGIN_FileWriter__init__

#------------------------------------------------------------------------------------

# Replace writers.CPPWriter by PLUGIN_FileWriter (remove formatting)
class PLUGIN_FileWriter(PLUGIN_writers.FileWriter):
    """Default FileWriter with minimal modifications"""

DEFAULT_CPPWriter = PLUGIN_writers.CPPWriter
#PLUGIN_writers.CPPWriter = DEFAULT_FileWriter # WITH FORMATTING
PLUGIN_writers.CPPWriter = PLUGIN_FileWriter # WITHOUT FORMATTING

#------------------------------------------------------------------------------------

import aloha
import aloha.aloha_writers as aloha_writers

from collections import defaultdict
from fractions import Fraction
from six import StringIO

# Define a custom ALOHAWriter
# (NB: enable this via PLUGIN_UFOModelConverter.aloha_writer)
class PLUGIN_ALOHAWriter(aloha_writers.ALOHAWriterForGPU):
    # Class structure information
    #  - object
    #  - WriteALOHA(object) [in aloha/aloha_writers.py]
    #  - ALOHAWriterForCPP(WriteALOHA) [in aloha/aloha_writers.py]
    #  - ALOHAWriterForGPU(ALOHAWriterForCPP) [in aloha/aloha_writers.py]
    #  - PLUGIN_ALOHAWriter(ALOHAWriterForGPU)
    #      This class

    # Keep defaults from aloha_writers.ALOHAWriterForGPU
    #extension = '.cc'
    prefix ='SYCL_EXTERNAL'
    type2def = {}
    type2def['pointer_vertex'] = '*' # using complex<double>* vertex
    type2def['pointer_coup'] = ''

    # Modify C++ code from aloha_writers.ALOHAWriterForGPU
    #ci_definition = 'cxtype cI = cxtype(0., 1.);\n'
    #ci_definition = 'const cxtype cI = cxmake( 0., 1. );\n'
    #realoperator = '.real()'
    #imagoperator = '.imag()'
    realoperator = 'CXREAL' # N.B. now a preprocessor macro 
    imagoperator = 'CXIMAG' # N.B. now a preprocessor macro 

    # Improve formatting
    type2def['int'] = 'int'
    type2def['double'] = 'fptype'
    type2def['complex'] = 'cxtype'

    # Add vector types
    type2def['double_v'] = 'fptype_sv'
    type2def['complex_v'] = 'cxtype_sv'

    # Modify C++ code from aloha_writers.ALOHAWriterForGPU
    # New option: declare C++ variable type only when they are defined?
    #nodeclare = False # old behaviour (separate declaration with no initialization)
    nodeclare = True # new behaviour (delayed declaration with initialisation)

    # Modify aloha_writers.ALOHAWriterForCPP method (improve formatting)
    def change_number_format(self, number):
        """Formatting the number"""
        def isinteger(x):
            try:
                return int(x) == x
            except TypeError:
                return False
        if isinteger(number):
            if number == 1: out = 'one'
            elif number == -1: out = '- one'
            elif number == 2: out = 'two'
            elif number == -2: out = '- two'
            else: out = '%s.' % (str(int(number))) # This prints -1 as '-1.'
        elif isinstance(number, complex):
            if number.imag:
                if number.real:
                    out = '( %s + %s * CXIMAGINARYI_SV )' % (self.change_number_format(number.real), self.change_number_format(number.imag))
                else:
                    if number.imag == 1:
                        out = 'CXIMAGINARYI_SV'
                    elif number.imag == -1:
                        out = '- CXIMAGINARYI_SV'
                    else:
                        out = '%s * CXIMAGINARYI_SV' % self.change_number_format(number.imag)
            else:
                out = '%s' % (self.change_number_format(number.real))
        else:
            tmp = Fraction(str(number))
            tmp = tmp.limit_denominator(100)
            if not abs(tmp - number) / abs(tmp + number) < 1e-8: out = '%.9f' % (number)
            elif tmp.numerator == 1 and tmp.denominator == 2 : out = 'half' # AV
            elif tmp.numerator == -1 and tmp.denominator == 2 : out = '- half' # AV
            else: out = '%s./%s.' % (tmp.numerator, tmp.denominator)
        return out

    # Modify aloha_writers.ALOHAWriterForCPP method (improve formatting)
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
                #type = self.type2def[format[5:]] # double or complex (instead of list_double or list_complex)
                #if not argname.startswith('COUP'): type = self.type2def[format[5:]+'_v'] # vectorize (double_v or complex_v)
                type = self.type2def[format[5:]+'_v'] # vectorize (double_v or complex_v)
                list_arg = '[]'
                comment_inputs.append('%s[6]'%argname) # N.B. (wavefuncsize=6 is hardcoded also in export_cpp...)
            else:
                type = self.type2def[format]
                list_arg = ''
            if argname.startswith('COUP'):
                type += '_sv'
                point = self.type2def['pointer_coup']
                args.append('%s %s%s%s'% (type, point, argname, list_arg))
            else:
                args.append('%s %s%s'% (type, argname, list_arg))
        if not self.offshell:
            #output = '%(doublec)s %(pointer_vertex)s vertex' % {
            output = '%(doublec)s%(pointer_vertex)s vertex' % { # vectorize
                'doublec':self.type2def['complex_v'],
                'pointer_vertex': self.type2def['pointer_vertex']}
            comment_output = 'amplitude \'vertex\''
        else:
            output = '%(doublec)s %(spin)s%(id)d[]' % { # vectorize
                     'doublec': self.type2def['complex_v'],
                     'spin': self.particles[self.outgoing -1],
                     'id': self.outgoing}
            comment_output = 'wavefunction \'%s%d[6]\'' % ( self.particles[self.outgoing -1], self.outgoing ) # (wavefuncsize=6)
        comment = '// Compute the output %s from the input wavefunctions %s' % ( comment_output, ', '.join(comment_inputs) )
        indent = ' ' * len( '  void %s( ' % name )
        out.write('  %(comment)s\n  %(prefix)s\n  void %(name)s( const %(args)s,\n%(indent)s%(output)s )%(suffix)s' %
                  {'comment': comment, # AV - add comment
                   'prefix': self.prefix + ( ' INLINE' if 'is_h' in mode else '' ), # AV - add INLINE
                   'suffix': ( ' ALWAYS_INLINE' if 'is_h' in mode else '' ), # AV - add ALWAYS_INLINE
                   'indent':indent, 'output':output, 'name': name,
                   'args': (',\n' + indent + 'const ').join(args)}) # AV - add const, add indent
        if 'is_h' in mode:
            out.write(';\n')
            out.write('\n  //--------------------------------------------------------------------------\n') # add footer
        else:
            out.write('\n  {\n')
        return out.getvalue()

    # Modify aloha_writers.ALOHAWriterForCPP method (improve formatting)
    # This affects HelAmps_sm.cc
    def get_foot_txt(self):
        """Prototype for language specific footer"""
        return '  }\n\n  //--------------------------------------------------------------------------'

    # Modify aloha_writers.ALOHAWriterForCPP method (improve formatting)
    # This affects HelAmps_sm.cc
    def get_declaration_txt(self, add_i=True):
        """ Prototype for how to write the declaration of variable
            Include the symmetry line (entry FFV_2)
        """
        out = StringIO()
        argument_var = [name for type,name in self.call_arg]
        # define the complex number CI = 0+1j
        if add_i:
            #out.write(self.ci_definition)
            #out.write('    ' + self.ci_definition)
            out.write('')
        codedict = {} # Allow delayed declaration with initialisation
        for type, name in self.declaration.tolist():
            #print(name) # FOR DEBUGGING
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
                fullname = '%s[%s]'%(name, size)
            elif (type, name) not in self.call_arg:
                fullname = name
            else:
                continue # No need to declare the variable
            if fullname.startswith('OM') :
                codedict[fullname] = '%s %s' % (self.type2def[type], fullname) # AV UGLY HACK (OM3 is always a scalar)
            else:
                codedict[fullname] = '%s %s' % (self.type2def[type+"_v"], fullname) # AV vectorize, add to codedict
            ###print(fullname, codedict[fullname]) # FOR DEBUGGING
            if self.nodeclare:
                self.declaration.codedict = codedict # New behaviour (delayed declaration with initialisation)
            else:
                out.write('    %s;\n' % codedict[fullname] ) # Old behaviour (separate declaration with no initialization)
        ###out.write('    // END DECLARATION\n') # FOR DEBUGGING
        return out.getvalue()

    # Modify aloha_writers.ALOHAWriterForCPP method (improve formatting)
    # This affects 'V1[0] = ' in HelAmps_sm.cc
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
                #out.write("    OM{0} = {1};\n    if (M{0} != {1})\n OM{0}={2}/(M{0}*M{0});\n".format(
                #out.write("    OM{0} = {1};\n    if ( M{0} != {1} ) OM{0} = {2} / (M{0}*M{0});\n".format(
                #out.write("    OM{0} = ( M{0} != {1} ? {2} / ( M{0} * M{0} ) : {1} );\n".format( # Use ternary in OM3
                #    #i+1, self.change_number_format(0), self.change_number_format(1)))
                #    i+1, '0.', '1.')) # AV force scalar "1." instead of vector "one"
                declname = 'OM%s' % (i+1) # AV
                if self.nodeclare: declname = 'const ' + self.declaration.codedict[declname]
                out.write("    {3} = ( M{0} != {1} ? {2} / ( M{0} * M{0} ) : {1} );\n".format( # Use ternary in OM3
                    i+1, '0.', '1.', declname)) # Force scalar "1." instead of vector "one", add declaration
            if i+1 == self.outgoing:
                out_type = type
                out_size = self.type_to_size[type]
                continue
            elif self.offshell:
                if len(p) != 0 : p.append(' ')
                #p.append('{0}{1}{2}[%(i)s]'.format(signs[i],type,i+1,type))
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
                #out.write('    %s%s[%s] = %s;\n' % (type,self.outgoing, i, ''.join(p) % dict_energy))
                out.write('    %s%s[%s] = %s;\n' % (type,self.outgoing, i, ''.join(p) % dict_energy))
            if self.declaration.is_used('P%s' % self.outgoing):
                self.get_one_momenta_def(self.outgoing, out)
        # Returning result
        return out.getvalue()

    # Modify aloha_writers.ALOHAWriterForCPP method (improve formatting, add delayed declaration with initialisation)
    # This affects 'P1[0] = ' in HelAmps_sm.cc
    def get_one_momenta_def(self, i, strfile):
        type = self.particles[i-1]
        if aloha.loop_mode:
            ptype = 'complex_v'
            templateval ='%(sign)s %(type)s%(i)d[%(nb)d]'
        else:
            ptype = 'double_v'
            templateval ='%(sign)s%(operator)s( %(type)s%(i)d[%(nb2)d] )' # cxreal/cximag
        if self.nodeclare: strfile.write('    const %s P%d[4] = { ' % ( self.type2def[ptype], i) )
        nb2 = 0
        for j in range(4):
            if not aloha.loop_mode:
                nb = j
                if j == 0:
                    misc.sprint('assert not aloha.mp_precision')
                    assert not aloha.mp_precision
                    operator = self.realoperator # not suppose to pass here in mp
                elif j == 1:
                    nb2 += 1
                elif j == 2:
                    misc.sprint('assert not aloha.mp_precision')
                    assert not aloha.mp_precision
                    operator = self.imagoperator # not suppose to pass here in mp
                elif j ==3:
                    nb2 -= 1
            else:
                operator =''
                nb = j
                nb2 = j
            sign = self.get_P_sign(i) if self.get_P_sign(i) else '+'
            if self.nodeclare: template = templateval + ( ', ' if j<3 else '' )
            else: template ='    P%(i)d[%(j)d] = ' + templateval + ';\n'
            strfile.write(template % {'j':j,'type': type, 'i': i,
                        'nb': nb, 'nb2': nb2, 'operator':operator,
                        'sign': sign})
        if self.nodeclare: strfile.write(' };\n')

    # Modify aloha_writers.ALOHAWriterForCPP method (improve formatting)
    # This is called once per FFV function, i.e. once per WriteALOHA instance?
    # It is called by WriteALOHA.write, after get_header_txt, get_declaration_txt, get_momenta_txt, before get_foot_txt
    # This affects 'denom = COUP' in HelAmps_sm.cc
    # This affects 'V1[2] = ' and 'F1[2] = ' in HelAmps_sm.cc
    # This affects 'TMP0 = ' in HelAmps_sm.cc
    # This affects '(*vertex) = ' in HelAmps_sm.cc
    def define_expression(self):
        """Write the helicity amplitude in C++ format"""
        out = StringIO()
        if self.routine.contracted:
            keys = sorted(self.routine.contracted.keys())
            for name in keys:
                obj = self.routine.contracted[name]
                # This affects 'TMP0 = ' in HelAmps_sm.cc
                if self.nodeclare:
                    out.write('    const %s %s = %s;\n' %
                              (self.type2def['complex_v'], name, self.write_obj(obj)))
                else:
                    out.write('    %s = %s;\n' % (name, self.write_obj(obj)))
                    self.declaration.add(('complex', name))
        for name, (fct, objs) in self.routine.fct.items():
            format = ' %s = %s;\n' % (name, self.get_fct_format(fct))
            out.write(format % ','.join([self.write_obj(obj) for obj in objs])) # N.B. not used in eemumu?
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
                # This affects '(*vertex) = ' in HelAmps_sm.cc
                out.write('    %(pre_vertex)svertex%(post_vertex)s = %(pre_coup)sCOUP%(post_coup)s * %(num)s;\n' % mydict)
            else:
                mydict= {}
                if self.type2def['pointer_vertex'] in ['*']:
                    mydict['pre_vertex'] = '(*'
                    mydict['post_vertex'] = ')'
                else:
                    mydict['pre_vertex'] = ''
                    mydict['post_vertex'] = ''
                mydict['data'] = self.write_obj(numerator.get_rep([0]))
                # This affects '(*vertex) = ' in HelAmps_sm.cc
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
                    mydict['declnamedenom'] = 'const %s denom' % self.type2def['complex_v']
                else:
                    mydict['declnamedenom'] = 'denom'
                    self.declaration.add(('complex','denom'))
                if not aloha.complex_mass:
                    # This affects 'denom = COUP' in HelAmps_sm.cc
                    if self.routine.denominator:
                        out.write('    %(declnamedenom)s = %(pre_coup)s%(coup)s%(post_coup)s / (%(denom)s)\n' % mydict)
                    else:
                        out.write('    %(declnamedenom)s = %(pre_coup)s%(coup)s%(post_coup)s / ( (P%(i)s[0] * P%(i)s[0] ) - ( P%(i)s[1] * P%(i)s[1] ) - ( P%(i)s[2] * P%(i)s[2] ) - ( P%(i)s[3] * P%(i)s[3] ) - M%(i)s * ( M%(i)s - CXIMAGINARYI_SV * W%(i)s ) );\n' % mydict)
                else:
                    if self.routine.denominator:
                        raise Exception('modify denominator are not compatible with complex mass scheme')
                    # This affects 'denom = COUP' in HelAmps_sm.cc
                    out.write('    %(declnamedenom)s = %(pre_coup)s%(coup)s%(post_coup)s / ( (P%(i)s[0] * P%(i)s[0] ) - ( P%(i)s[1] *P%(i)s[1] ) - ( P%(i)s[2] * P%(i)s[2] ) - ( P%(i)s[3] * P%(i)s[3] ) - ( M%(i)s * M%(i)s ) );\n' % mydict)
                ###self.declaration.add(('complex','denom')) # AV moved earlier (or simply removed)
                if aloha.loop_mode: ptype = 'list_complex'
                else: ptype = 'list_double'
                self.declaration.add((ptype,'P%s' % self.outgoing))
            else:
                coeff = 'COUP'
            for ind in numerator.listindices():
                # This affects 'V1[2] = ' and 'F1[2] = ' in HelAmps_sm.cc
                out.write('    %s[%d] = %s * %s;\n' % (self.outname, self.pass_to_HELAS(ind), coeff, self.write_obj(numerator.get_rep(ind))))
        # FIXME check if one, two or half are used and need to be defined (ugly hack for #291: can this be done better?)
        out2 = StringIO()
        if 'one' in out.getvalue(): out2.write('    constexpr fptype one( 1. );\n')
        if 'two' in out.getvalue(): out2.write('    constexpr fptype two( 2. );\n')
        if 'half' in out.getvalue(): out2.write('    constexpr fptype half( 1. / 2. );\n')
        out2.write( out.getvalue() )
        return out2.getvalue()

    # Modify aloha_writers.WriteALOHA method (improve formatting)
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

    # Modify aloha_writers.WriteALOHA method (improve formatting)
    def write_MultContainer(self, obj, prefactor=True):
        """Turn a multvariable into a string"""
        mult_list = [self.write_obj(id) for id in obj]
        data = {'factors': ' * '.join(mult_list)}
        if prefactor and obj.prefactor != 1:
            if obj.prefactor != -1:
                text = '%(prefactor)s * %(factors)s'
                data['prefactor'] = self.change_number_format(obj.prefactor)
            else:
                text = '-%(factors)s' # Keep default (this is not used in eemumu)
        else:
            text = '%(factors)s'
        return text % data

    # Modify aloha_writers.WriteALOHA method (improve formatting)
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
                file_str.write('(%s)' % formatted)
            else:
                file_str.write(formatted)
            file_str.write(' * ( ')
        else:
            file_str.write('( ')
        first=True
        for value, obj_list in data.items():
            add= ' + '
            if value not in  [-1,1]:
                nb_str = self.change_number_format(value)
                if nb_str[0] in ['+','-']:
                    file_str.write(nb_str)
                else:
                    file_str.write('+' if first else ' + ')
                    file_str.write(nb_str)
                file_str.write(' * ( ') # (eg '+ CXIMAGINARYI_SV * (V3[4])')
            elif value == -1:
                add = ' - '
                file_str.write('-' if first else ' - ')
            elif not first:
                file_str.write(' + ')
            else:
                file_str.write('')
            first = False
            # N.B. write_obj here also adds calls declaration_add (via change_var_format) - example: OM3
            file_str.write(add.join([self.write_obj(obj, prefactor=False) for obj in obj_list]))
            if value not in [1,-1]:
                file_str.write(' )')
        if number:
            total = sum(number)
            file_str.write('+ %s' % self.change_number_format(total))
        file_str.write(' )')
        return file_str.getvalue()

#------------------------------------------------------------------------------------

class PLUGIN_UFOModelConverter(PLUGIN_export_cpp.UFOModelConverterGPU):
    # Class structure information
    #  - object
    #  - UFOModelConverterCPP(object) [in madgraph/iolibs/export_cpp.py]
    #  - UFOModelConverterGPU(UFOModelConverterCPP) [in madgraph/iolibs/export_cpp.py]
    #  - PLUGIN_UFOModelConverter(UFOModelConverterGPU)
    #      This class

    # Keep defaults from export_cpp.UFOModelConverterCPP
    #include_dir = '.'
    #c_file_dir = '.'
    #param_template_h = 'cpp_model_parameters_h.inc'
    #param_template_cc = 'cpp_model_parameters_cc.inc'

    # Change defaults from export_cpp.UFOModelConverterCPP
    # (custom tag to appear in 'This file has been automatically generated for')
    output_name = 'SYCL standalone'

    # Change defaults from export_cpp.UFOModelConverterGPU
    cc_ext = 'cc' # create HelAmps_sm.cc

    # Keep defaults from export_cpp.UFOModelConverterGPU
    #cc_ext = 'cc'
    #aloha_template_h = pjoin('gpu','cpp_hel_amps_h.inc')
    #aloha_template_cc = pjoin('gpu','cpp_hel_amps_cc.inc')
    #helas_h = pjoin('gpu', 'helas.h')
    #helas_cc = pjoin('gpu', 'helas.cu')

    # Use a custom ALOHAWriter (NB: this is an argument to WriterFactory.__new__, either a string or a class!)
    aloha_writer = PLUGIN_ALOHAWriter # WriterFactory will use ALOHAWriterForGPU

    # Use template files from PLUGINDIR instead of MG5DIR
    def read_aloha_template_files(self, ext):
        """Read all ALOHA template files with extension ext, strip them of
        compiler options and namespace options, and return in a list"""
        path = pjoin(PLUGINDIR, 'aloha', 'template_files')
        out = []
        file = ''
        if ext == 'h':
            file = open(pjoin(path, self.helas_h)).read()
            file = '\n'.join( file.split('\n')[12:] ) # skip first 12 lines in helas.h (copyright)
        else:
            file = open(pjoin(path, self.helas_cc)).read()
            file = '\n'.join( file.split('\n')[8:] ) # skip first 8 lines in helas.cu (copyright)
        out.append( file )
        return out

    # Use the plugin's PLUGIN_OneProcessExporter template_path and __template_path (for aloha_template_h/cc)
    @classmethod
    def read_template_file(cls, filename, classpath=False):
        """Open a template file and return the contents."""
        return PLUGIN_OneProcessExporter.read_template_file(filename, classpath)

    # Overload export_cpp.UFOModelConverterCPP method (improve formatting)
    def write_parameters(self, params):
        res = super().write_parameters(params)
        res = res.replace('std::complex<','mgOnGpu::cxsmpl<') # custom simplex complex class (with constexpr arithmetics)
        if res == '' : res = '  // (none)'
        else : res = '  ' + res # add leading '  ' after the '// Model' line
        res = res.replace('\n','\n  ')
        res = res.replace(',',', ')
        return res

    # Overload export_cpp.UFOModelConverterCPP method (improve formatting)
    def write_set_parameters(self, params):
        res = self.super_write_set_parameters_donotfixMajorana(params)
        res = res.replace('std::complex<','mgOnGpu::cxsmpl<') # custom simplex complex class (with constexpr arithmetics)
        if res == '' : res = '// (none)'
        res = res.replace('\n','\n  ')
        return res

    def write_hardcoded_parameters(self, params):
        pardef = super().write_parameters(params)
        parset = self.super_write_set_parameters_donotfixMajorana(params)
        if ( pardef == '' ):
            misc.sprint('assert parset == \'\'')
            assert parset == '', "pardef is empty but parset is not: '%s'"%parset # AV sanity check (both are empty)
            res = '// (none)\n'
            return res
        pardef = pardef.replace('std::complex<','mgOnGpu::cxsmpl<') # custom simplex complex class (with constexpr arithmetics)
        parset = parset.replace('std::complex<','mgOnGpu::cxsmpl<') # custom simplex complex class (with constexpr arithmetics)
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
        misc.sprint( 'pardef_lines size =', len(pardef_lines), ', keys size =', len(pardef_lines.keys()) )

        parset_pars = []
        parset_lines = {}
        skipnextline = False
        for iline, line in enumerate(parset.split('\n')):
        for line in parset.split('\n'):
            if line.startswith('indices'):
                continue # skip line with leading "indices", before slha.get_block_entry (#622)
            par, parval = line.split(' = ')
            if parval.startswith('slha.get_block_entry'): parval = parval.split(',')[2].lstrip(' ').rstrip(');') + ';'
            parset_pars.append( par )
            parset_lines[par] = parval # includes a trailing ';'
        misc.sprint( 'parset_pars size =', len(parset_pars) )
        misc.sprint( 'parset_lines size =', len(parset_lines), ', keys size =', len(parset_lines.keys()) )
        assert len(pardef_lines) == len(parset_lines), 'len(pardef_lines) != len(parset_lines)' # sanity check (same number of parameters)

        res = '  '.join( pardef_lines[par] + ' = ' + parset_lines[par] + '\n' for par in parset_pars ) # no leading '  ' on first row
        res = res.replace(' ;',';')
        res = res.replace('= - ','= -') # post-fix for susy
        res = res.replace('(  - ','( -') # post-fix for susy
        res = res.replace(',  - ',', -') # post-fix for SM no_b_mass
        return res

    # replace export_cpp.UFOModelConverterCPP method (split writing of parameters and fixes for Majorana particles #622)
    def super_write_set_parameters_donotfixMajorana(self, params):
        """Write out the lines of independent parameters"""
        res_strings = []
        # For each parameter, write name = expr;
        for param in params:
            res_strings.append("%s" % param.expr)
        return "\n".join(res_strings)

    # replace export_cpp.UFOModelConverterCPP method (eventually split writing of parameters and fixes for Majorana particles #622)
    def super_write_set_parameters_onlyfixMajorana(self, hardcoded): # FIXME! split hardcoded (constexpr) and not-hardcoded code
        """Write out the lines of independent parameters"""
        print( 'super_write_set_parameters_onlyfixMajorana (hardcoded=%s)'%hardcoded )
        res_strings = []
        # Correct width sign for Majorana particles (where the width and mass need to have the same sign)        
        for particle in self.model.get('particles'):
            if particle.is_fermion() and particle.get('self_antipart') and \
                   particle.get('width').lower() != 'zero':
                res_strings.append("  if( %s < 0 )" % particle.get('mass'))
                res_strings.append("    %(width)s = -abs( %(width)s );" % {"width": particle.get('width')})
        return '\n' + '\n'.join(res_strings) if res_strings else ''


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
        replace_dict['print_independent_parameters'] += self.super_write_set_parameters_onlyfixMajorana( hardcoded=False ) # add fixes for Majorana particles only in the aS-indep parameters #622
        replace_dict['print_independent_couplings'] = self.write_print_parameters(self.coups_indep)
        replace_dict['print_dependent_parameters'] = self.write_print_parameters(self.params_dep)
        replace_dict['print_dependent_couplings'] = self.write_print_parameters(list(self.coups_dep.values()))
        if 'include_prefix' not in replace_dict:
            replace_dict['include_prefix'] = ''
        assert super().write_parameters([]) == '', 'super().write_parameters([]) is not empty' # sanity check (#622)
        assert self.super_write_set_parameters_donotfixMajorana([]) == '', 'super_write_set_parameters_donotfixMajorana([]) is not empty' # sanity check (#622)
        hrd_params_indep = [ line.replace('constexpr','//constexpr') + ' // now retrieved event-by-event (as G) from Fortran (running alphas #373)' if 'aS =' in line else line for line in self.write_hardcoded_parameters(self.params_indep).split('\n') ]
        replace_dict['hardcoded_independent_parameters'] = '\n'.join( hrd_params_indep ) + self.super_write_set_parameters_onlyfixMajorana( hardcoded=True ) # add fixes for Majorana particles only in the aS-indep parameters #622
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
            dcoupdecl = [ '    cxtype_sv %s;' % name for name in self.coups_dep ]
            replace_dict['dcoupdecl'] = '\n'.join( dcoupdecl )
            dcoupsetdpar = []
            foundG = False
            for line in self.write_hardcoded_parameters(self.params_dep).split('\n'):
                if line != '':
                    dcoupsetdpar.append( '    ' + line.replace('constexpr double', 'const FPType' if foundG else '//const FPType' ) )
                    if 'constexpr double G =' in line: foundG = True
            replace_dict['dcoupsetdpar'] = '  ' + '\n'.join( dcoupsetdpar )
            dcoupsetdcoup = [ '    ' + line.replace('constexpr mgOnGpu::cxsmpl<double> ','couplings[idcoup_')
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

    # Overload export_cpp.UFOModelConverterCPP method (improve formatting)
    def generate_parameters_class_files(self):
        file_h, file_cc = self.super_generate_parameters_class_files()
        file_h = file_h[:-1] # remove extra trailing '\n'
        file_cc = file_cc[:-1] # remove extra trailing '\n'
        # [NB: there is a minor bug in export_cpp.UFOModelConverterCPP.generate_parameters_class_files
        # ['independent_couplings' contains dependent parameters, 'dependent parameters' contains independent_couplings]
        # [This only affects the order in which they are printed out - which is now reversed in the templates]
        return file_h, file_cc

    # Replace export_cpp.UFOModelConverterCPP method (add explicit std namespace)
    def write_print_parameters(self, params):
        """Write out the lines of independent parameters"""
        # For each parameter, write name = expr;
        res_strings = []
        for param in params:
            res_strings.append('std::cout << std::setw( 20 ) << \"%s = \" << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << %s << std::endl;' % (param.name, param.name)) # AV
        if len(res_strings) == 0 : res_strings.append('// (none)')
        return "\n  ".join(res_strings) # FIXME (why was this not necessary before?)

    # Replace export_cpp.UFOModelConverterCPP method (add debug printouts)
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
        """Generate the hel_amps_model.h and hel_amps_model.cc files that have the complete set of generalized Helas routines for the model"""
        import aloha.create_aloha as create_aloha
        if not os.path.isdir(os.path.join(self.dir_path, self.include_dir)):
            os.makedirs(os.path.join(self.dir_path, self.include_dir))
        if not os.path.isdir(os.path.join(self.dir_path, self.cc_file_dir)):
            os.makedirs(os.path.join(self.dir_path, self.cc_file_dir))
        model_h_file  = os.path.join(self.dir_path, self.include_dir, 'HelAmps_%s.h'  % self.model_name)
        model_cc_file = os.path.join(self.dir_path, self.cc_file_dir, 'HelAmps_%s.%s' % (self.model_name, self.cc_ext))
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
            print(type(abstracthelas), abstracthelas.name) # N.B. this is the loop on FFV functions
            h_rout, cc_rout = abstracthelas.write(output_dir=None, language=self.aloha_writer, mode='no_include')
            template_h_files.append(h_rout)
            template_cc_files.append(cc_rout)
        replace_dict['function_declarations'] = '\n'.join(template_h_files)
        replace_dict['function_definitions'] = '\n'.join(template_cc_files)
        file_h = self.read_template_file(self.aloha_template_h) % replace_dict
        file_cc = self.read_template_file(self.aloha_template_cc) % replace_dict
        file_cc = '\n'.join( file_cc.split('\n')[8:] ) # skip first 8 lines in cpp_hel_amps_cc.inc (copyright)
        # Write the files
        file_h_lines = file_h.split('\n')
        file_h = '\n'.join( file_h_lines[:-3]) # skip the trailing '//---'
        file_h += file_cc # append the contents of HelAmps_sm.cc directly to HelAmps_sm.h!
        PLUGIN_writers.CPPWriter(model_h_file).writelines(file_h)
        #PLUGIN_writers.CPPWriter(model_cc_file).writelines(file_cc)
        #logger.info("Created files %s and %s in directory" \
        #            % (os.path.split(model_h_file)[-1],
        #               os.path.split(model_cc_file)[-1]))
        logger.info("Created file %s in directory %s" \
                    % (os.path.split(model_h_file)[-1], os.path.split(model_h_file)[0] ) )

#------------------------------------------------------------------------------------

import madgraph.iolibs.files as files
import madgraph.various.misc as misc

class PLUGIN_OneProcessExporter(PLUGIN_export_cpp.OneProcessExporterGPU):
    # Class structure information
    #  - object
    #  - OneProcessExporterCPP(object) [in madgraph/iolibs/export_cpp.py]
    #  - OneProcessExporterGPU(OneProcessExporterCPP) [in madgraph/iolibs/export_cpp.py]
    #  - PLUGIN_OneProcessExporter(OneProcessExporterGPU)
    #      This class

    # Change defaults from export_cpp.OneProcessExporterGPU
    # [NB process_class = "CPPProcess" is set in OneProcessExporterCPP.__init__]
    cc_ext = 'cc' # create CPPProcess.cc

    # Keep defaults from export_cpp.OneProcessExporterGPU
    #process_dir = '.'
    #include_dir = '.'
    #process_template_h = 'gpu/process_h.inc'
    #process_template_cc = 'gpu/process_cc.inc'
    #process_class_template = 'gpu/process_class.inc'
    #process_definition_template = 'gpu/process_function_definitions.inc'
    #process_wavefunction_template = 'cpp_process_wavefunctions.inc'
    #process_sigmaKin_function_template = 'gpu/process_sigmaKin_function.inc'
    #single_process_template = 'gpu/process_matrix.inc'

    # Use template files from MG5DIR
    #template_path = os.path.join(_file_path, 'iolibs', 'template_files')
    #__template_path = os.path.join(_file_path, 'iolibs', 'template_files')

    # Use template files from PLUGINDIR instead of MG5DIR
    template_path = os.path.join( PLUGINDIR, 'madgraph', 'iolibs', 'template_files' )
    __template_path = os.path.join( PLUGINDIR, 'madgraph', 'iolibs', 'template_files' )

    # Overload export_cpp.OneProcessExporterGPU constructor (rename gCPPProcess to CPPProcess, set include_multi_channel)
    def __init__(self, *args, **kwargs):
        misc.sprint('Entering PLUGIN_OneProcessExporter.__init__')
        for kwarg in kwargs: misc.sprint( 'kwargs[%s] = %s' %( kwarg, kwargs[kwarg] ) )
        super().__init__(*args, **kwargs)
        self.process_class = 'CPPProcess'
        if 'prefix' in kwargs: proc_id = kwargs['prefix']+1 # madevent+sycl (ime+1 from ProcessExporterFortranMEGroup.generate_subprocess_directory)
        else: proc_id = 0 # standalone_sycl
        misc.sprint(proc_id)
        self.proc_id = proc_id

    # Modify export_cpp.OneProcessExporterGPU method (indent comments in process_lines)
    def get_process_class_definitions(self, write=True):
        replace_dict = super().get_process_class_definitions(write=False)
        replace_dict['process_lines'] = replace_dict['process_lines'].replace('\n','\n  ')
        replace_dict['nwavefunc'] = self.matrix_elements[0].get_number_of_wavefunctions()
        replace_dict['wavefuncsize'] = 6
        nexternal, nincoming = self.matrix_elements[0].get_nexternal_ninitial()
        replace_dict['nincoming'] = nincoming
        replace_dict['noutcoming'] = nexternal - nincoming
        replace_dict['nexternal'] = nexternal
        replace_dict['nbhel'] = self.matrix_elements[0].get_helicity_combinations() # number of helicity combinations
        file = self.read_template_file(self.process_class_template) % replace_dict # HACK! ignore write=False case
        file = '\n'.join( file.split('\n')[12:] ) # skip first 12 lines in process_class.inc (copyright)
        return file

    # Modify export_cpp.OneProcessExporterGPU method (fix CPPProcess.cc)
    def get_process_function_definitions(self, write=True):
        misc.sprint('Entering PLUGIN_OneProcessExporter.get_process_function_definitions')
        """The complete class definition for the process"""
        #FIXME THINGS BREAK HERE FOR SUSY MODEL
        replace_dict = super(PLUGIN_export_cpp.OneProcessExporterGPU,self).get_process_function_definitions(write=False)
        replace_dict['ncouplings'] = len(self.couplings2order)
        replace_dict['ncouplingstimes2'] = 2 * replace_dict['ncouplings']
        replace_dict['nparams'] = len(self.params2order)
        replace_dict['nmodels'] = replace_dict['nparams'] + replace_dict['ncouplings']
        replace_dict['coupling_list'] = ' '
        replace_dict['hel_amps_cc'] = "#include \"HelAmps_%s.cc\"" % self.model_name # AV
        coupling = [''] * len(self.couplings2order)
        params = [''] * len(self.params2order)
        for coup, pos in self.couplings2order.items():
            coupling[pos] = coup

        coupling_indep = [] # Keep only the alphas-independent couplings #434
        for coup in coupling:
            keep = True
            # Use the same implementation as in UFOModelConverterCPP.prepare_couplings (assume self.model is the same)
            for key, coup_list in self.model['couplings'].items():
                if "aS" in key and coup in coup_list: keep = False
            if keep: coupling_indep.append( coup ) # N.B. only indep!

        # Need access to independent couplings tIPC outside of CPPProcess for SYCL
        coup_str = ""
        if len(coupling_indep) > 0:
            for i in range(len(coupling_indep)):
                coup_str += "          m_tIPC[%s] = cxmake( m_pars->%s );\n" % (i, coupling_indep[i])
        else:
            coup_str = "          //m_tIPC[...] = ... ; // nicoup=0\n"

        self.number_dependent_couplings = len(coupling) - len(coupling_indep)
        self.number_independent_couplings = len(coupling_indep)

        for para, pos in self.params2order.items():
            params[pos] = para

        # Need access to tIPD outside of CPPProcess for SYCL
        param_str = ""
        if len(params) > 0:
            for i in range(len(self.params2order)):
                param_str += "          m_tIPD[%s] = (fptype)m_pars->%s;\n" % (i, params[i])
        else:
            param_str += "          //m_tIPD[...] = ... ; // nparam=0\n"

        replace_dict['assign_coupling'] = coup_str + param_str
        color_amplitudes = [me.get_color_amplitudes() for me in self.matrix_elements] # as in OneProcessExporterCPP.get_process_function_definitions
        replace_dict['ncolor'] = len(color_amplitudes[0])
        file = self.read_template_file(self.process_definition_template) % replace_dict
        file = '\n'.join( file.split('\n')[12:] ) # skip first 12 lines in process_function_definitions.inc (copyright)
        return file

    # Modify export_cpp.OneProcessExporterGPU method (add debug printouts for multichannel #342)
    def get_sigmaKin_lines(self, color_amplitudes, write=True):
        misc.sprint('Entering PLUGIN_OneProcessExporter.get_sigmaKin_lines')
        misc.sprint(self.include_multi_channel)
        misc.sprint(self.support_multichannel)
        replace_dict = super().get_sigmaKin_lines(color_amplitudes, write=False)
        replace_dict['proc_id'] = self.proc_id if self.proc_id>0 else 1
        replace_dict['proc_id_source'] = 'madevent + sycl exporter' if self.proc_id>0 else 'standalone_sycl'
        if write:
            file = self.read_template_file(self.process_sigmaKin_function_template) % replace_dict
            file = '\n'.join( file.split('\n')[12:] ) # skip first 12 lines in process_sigmaKin_function.inc (copyright)
            return file, replace_dict
        else:
            return replace_dict

    # Modify export_cpp.OneProcessExporterGPU method (fix CPPProcess.cc)
    def get_all_sigmaKin_lines(self, color_amplitudes, class_name):
        """Get sigmaKin_process for all subprocesses for gCPPProcess.cc"""
        ret_lines = []
        if self.single_helicities:
            ret_lines.append("""
  // Evaluate |M|^2 for each subprocess
  // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event(s)
  SYCL_EXTERNAL INLINE
  fptype_sv calculate_wavefunctions( const vector4* __restrict__ allmomenta,      // input: momenta as vector4 
                                     #ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                                         fptype_sv* __restrict__ allNumerators,   // output: multichannel numerators, running_sum_over_helicities
                                         fptype_sv* __restrict__ allDenominators, // output: multichannel denominators, running_sum_over_helicities
                                         const size_t channelId,                  // input: multichannel channel id (1 to #diagrams); 0 to disable channel enhancement
                                     #endif
                                     const signed char*  __restrict__ cHel,
                                     const cxtype_sv* __restrict__ COUPs,
                                     const fptype* __restrict__ cIPD,
                                     fptype_sv* jamp2_sv                          // output: jamp2[ncolor] for color choice (nullptr if disabled)
                                   ) {
      using namespace MG5_sm;
      fptype_sv allMEs = FPZERO_SV;
\n""")
            ret_lines.append("""
      // Local TEMPORARY variables for a subset of Feynman diagrams in the given SYCL event (ievt)
      // [NB these variables are reused several times (and re-initialised each time) within the same event or event page]
      cxtype_sv w_sv[CPPPROCESS_NWF][CPPPROCESS_NW6]; // particle wavefunctions within Feynman diagrams (CPPPROCESS_NW6 is often 6, the dimension of spin 1/2 or spin 1 particles)
      cxtype_sv amp_sv[1]; // invariant amplitude for one given Feynman diagram

      // Local variables for the given SYCL event (ievt)
      cxtype_sv jamp_sv[ncolor] = {}; // sum of the invariant amplitudes for all Feynman diagrams in the event or event page

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
        #to_add = [] # FIXME - what is this for? comment it out
        #to_add.extend([self.get_matrix_single_process(i, me, color_amplitudes[i], class_name) for i, me in enumerate(self.matrix_elements)])
        file_extend = []
        for i, me in enumerate(self.matrix_elements):
            file = self.get_matrix_single_process( i, me, color_amplitudes[i], class_name )
            file = '\n'.join( file.split('\n')[12:] ) # skip first 12 lines in process_matrix.inc (copyright)
            file_extend.append( file )
        ret_lines.extend( file_extend )
        return "\n".join(ret_lines)

    # Modify export_cpp.OneProcessExporterGPU method (replace '# Process' by '// Process')
    def get_process_info_lines(self, matrix_element):
        """Return info lines describing the processes for this matrix element"""
        return"\n".join([ "// " + process.nice_string().replace('\n', '\n// * ') \
                         for process in matrix_element.get('processes')])

    # Replace the export_cpp.OneProcessExporterGPU method (invert .cc/.cu, add debug printouts)
    def generate_process_files(self):
        """Generate mgOnGpuConfig.h, CPPProcess.cc, check_sa.cc, CPPProcess.h, check_sa.cc""" 
        misc.sprint('Entering PLUGIN_OneProcessExporter.generate_process_files')
        if self.include_multi_channel:
            misc.sprint('self.include_multi_channel is already defined: this is madevent+second_exporter mode')
        else:
            misc.sprint('self.include_multi_channel is not yet defined: this is standalone_sycl mode') # see issue #473
        if self.matrix_elements[0].get('has_mirror_process'):
            self.matrix_elements[0].set('has_mirror_process', False)
            self.nprocesses/=2

        #super(PLUGIN_export_cpp.OneProcessExporterGPU, self).generate_process_files()
        # N.B. Explicity copy method and generate cc file first (need params2order for h file)
        # Generate the .h and .cc files needed for C++, for the processes described by multi_matrix_element

        # Create the files
        if not os.path.isdir(os.path.join(self.path, self.process_dir)):
            os.makedirs(os.path.join(self.path, self.process_dir))
        filename = os.path.join(self.path, self.process_dir, '%s.%s' % (self.process_class, self.cc_ext))

        misc.sprint("Start writing process_cc files")
        self.write_process_cc_file(PLUGIN_writers.CPPWriter(filename))
        misc.sprint("Finished writing process_cc files")


        if not os.path.isdir(os.path.join(self.path, self.include_dir)):
            os.makedirs(os.path.join(self.path, self.include_dir))
        filename = os.path.join(self.path, self.include_dir, '%s.h' % self.process_class)

        self.write_process_h_file(PLUGIN_writers.CPPWriter(filename))


        logger.info('Created files %(process)s.h and %(process)s.cc in' % \
                    {'process': self.process_class} + \
                    ' directory %(dir)s' % {'dir': os.path.split(filename)[0]})

        self.edit_CMakeLists()
        self.edit_check_sa()
        self.edit_mgonGPU()
        self.edit_processidfile() # New file (N.B. this is Sigma-specific, should not be a symlink to Subprocesses)
        if self.include_multi_channel:
            self.edit_coloramps() # New file (N.B. this is Sigma-specific, should not be a symlink to Subprocesses)
        # Add symbolic links
        # N.B. symlink of sycl.mk to makefile is overwritten by madevent makefile if this exists (#480)
        # N.B. this relies on the assumption that sycl code is generated before madevent code
        files.ln(pjoin(self.path, 'sycl.mk'), self.path, 'makefile')

    # Generate CMakeLists.txt file inside the P* directory
    def edit_CMakeLists(self):
        """Generate CMakeLists.txt"""
        misc.sprint('Entering PLUGIN_OneProcessExporter.edit_CMakeLists')
        template = open(pjoin(self.template_path,'CMake/SubProcesses/CMakeLists_P.txt'),'r').read()
        ff = open(pjoin(self.path, 'CMakeLists.txt'),'w')
        ff.write(template)
        ff.close()

    # Replace the export_cpp.OneProcessExporterGPU method (invert .cc/.cu, add debug printouts)
    def edit_check_sa(self):
        """Generate check_sa.cc and fcheck_sa.f"""
        misc.sprint('Entering PLUGIN_OneProcessExporter.edit_check_sa')
        ff = open(pjoin(self.path, 'check_sa.cc'),'w')
        template = open(pjoin(self.template_path,'gpu','check_sa.cc'),'r').read()
        ff.write(template) # nothing to replace in check_sa.cc
        ff.close()
        replace_dict = {}
        replace_dict['nexternal'], _ = self.matrix_elements[0].get_nexternal_ninitial()
        #replace_dict['model'] = self.model_name
        #replace_dict['numproc'] = len(self.matrix_elements)
        ff = open(pjoin(self.path, 'fcheck_sa.f'),'w')
        template = open(pjoin(self.template_path,'gpu','fcheck_sa.f'),'r').read()
        ff.write(template % replace_dict)
        ff.close()

    # Add debug printouts over the export_cpp.OneProcessExporterGPU method
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
        replace_dict['nexternal'] = nexternal
        
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

    # New method
    def edit_processidfile(self):
        """Generate epoch_process_id.h"""
        misc.sprint('Entering PLUGIN_OneProcessExporter.edit_processidfile')
        template = open(pjoin(self.template_path,'gpu','epoch_process_id.h'),'r').read()
        replace_dict = {}
        replace_dict['processid'] = self.get_process_name().upper()
        replace_dict['processid_uppercase'] = self.get_process_name().upper()
        ff = open(pjoin(self.path, 'epoch_process_id.h'),'w')
        ff.write(template % replace_dict)
        ff.close()

    # New method
    def edit_coloramps(self):
        """Generate coloramps.h"""
        misc.sprint('Entering PLUGIN_OneProcessExporter.edit_coloramps')
        template = open(pjoin(self.template_path,'gpu','coloramps.h'),'r').read()
        ff = open(pjoin(self.path, 'coloramps.h'),'w')
        # The following five lines from OneProcessExporterCPP.get_sigmaKin_lines (using OneProcessExporterCPP.get_icolamp_lines)
        replace_dict={}
        if self.include_multi_channel: # N.B. unnecessary as edit_coloramps is not called otherwise...
            multi_channel = self.get_multi_channel_dictionary(self.matrix_elements[0].get('diagrams'), self.include_multi_channel)
            replace_dict['is_LC'] = self.get_icolamp_lines(multi_channel, self.matrix_elements[0], 1)
            #replace_dict['nb_channel'] = len(multi_channel)
            #replace_dict['nb_color'] = max(1,len(self.matrix_elements[0].get('color_basis')))
            # Extra formatting (e.g. gg_tt was "{{true,true};,{true,false};,{false,true};};")
            replace_dict['is_LC'] = replace_dict['is_LC'].replace(',',', ').replace('{{','      ').replace('};, {',',\n      ').replace('};};','')
        ff.write(template % replace_dict)
        ff.close()

    # Overload the export_cpp.OneProcessExporterGPU method (add debug printout and truncate last \n)
    # [*N.B. export_cpp.UFOModelConverterGPU.write_process_h_file is not called!*]
    def super_write_process_h_file(self, writer):
        """Generate final CPPProcess.h"""
        misc.sprint('Entering PLUGIN_OneProcessExporter.write_process_h_file')
        replace_dict = super(PLUGIN_export_cpp.OneProcessExporterGPU, self).write_process_h_file(False)
        #replace_dict2 = self.get_process_function_definitions(write=False)

        #Set helicities
        replace_dict['all_helicities'] = self.get_helicity_matrix(self.matrix_elements[0])
        replace_dict['all_helicities'] = replace_dict['all_helicities'] .replace("{", "")
        replace_dict['all_helicities'] = replace_dict['all_helicities'] .replace("}", "")
        replace_dict['all_helicities'] = replace_dict['all_helicities'] .replace("helicities", "constexpr T helicities[] {")
        replace_dict['all_helicities'] = replace_dict['all_helicities'] .replace("    ", "  ")
        replace_dict['all_helicities'] = replace_dict['all_helicities'] .replace(", constexpr", "constexpr")
        replace_dict['all_helicities'] = replace_dict['all_helicities'] .replace(";", "")

        #Set hardcoded parameters
        params = ['']
        params = [''] * len(self.params2order)
        for para, pos in self.params2order.items():
            params[pos] = para

        if len(params) > 0:
            param_str_hrd = '  constexpr FPType independent_parameters[] { '
            for para in params:
                param_str_hrd += '(FPType)Parameters_{0:s}::{1:s}, '.format( self.model_name, para )
            param_str_hrd = param_str_hrd[:-2] + ' };'
            replace_dict['independent_parameters_hrdcod'] = param_str_hrd
        else:
            replace_dict['independent_parameters_hrdcod'] = 'constexpr FPType independent_parameters[] {}; // unused as nparam=0'

        if writer:
            file = self.read_template_file(self.process_template_h) % replace_dict
            # Write the file
            writer.writelines(file)
        else:
            return replace_dict

    def write_process_h_file(self, writer):
        """Generate final CPPProcess.h"""
        misc.sprint('Entering PLUGIN_OneProcessExporter.write_process_h_file')
        #out = super().write_process_h_file(writer)
        out = self.super_write_process_h_file(writer)
        writer.seek(-1, os.SEEK_CUR)
        writer.truncate()
        return out

    # Replace the export_cpp.OneProcessExporterGPU method (replace HelAmps.cu by HelAmps.cc)
    def super_write_process_cc_file(self, writer):
        """Write the class member definition (.cc) file for the process described by matrix_element"""
        misc.sprint("Entering PLUGIN_OneProcessExporter.super_write_process_cc_file")
        replace_dict = super(PLUGIN_export_cpp.OneProcessExporterGPU, self).write_process_cc_file(False)
        replace_dict['hel_amps_h'] = "#include \"HelAmps_%s.h\"" % self.model_name
        if writer:
            file = self.read_template_file(self.process_template_cc) % replace_dict
            # Write the file
            writer.writelines(file)
        else:
            return replace_dict

    # Overload the export_cpp.OneProcessExporterGPU method (add debug printout and truncate last \n)
    def write_process_cc_file(self, writer):
        """Generate CPPProcess.cc"""
        misc.sprint('Entering PLUGIN_OneProcessExporter.write_process_cc_file')
        out = self.super_write_process_cc_file(writer)
        writer.seek(-1, os.SEEK_CUR)
        writer.truncate()
        return out

    # Replace the export_cpp.OneProcessExporterGPU method
    @staticmethod
    def coeff(ff_number, frac, is_imaginary, Nc_power, Nc_value=3):
        """Returns a nicely formatted string for the coefficients in JAMP lines"""
        total_coeff = ff_number * frac * Fraction(Nc_value) ** Nc_power
        if total_coeff == 1:
            if is_imaginary:
                return '+CXIMAGINARYI_SV*'
            else:
                return '+'
        elif total_coeff == -1:
            if is_imaginary:
                return '-CXIMAGINARYI_SV*'
            else:
                return '-' # Keep default (eg jamp_sv[0] += -amp_sv[0])
        res_str = '%+i.' % total_coeff.numerator
        if total_coeff.denominator != 1:
            # Check if total_coeff is an integer
            res_str = res_str + '/%i.' % total_coeff.denominator
        if is_imaginary:
            res_str = res_str + '*CXIMAGINARYI_SV'
        return res_str + '*'

    # Replace the export_cpp.OneProcessExporterCPP method (fix fptype and improve formatting)
    def get_color_matrix_lines(self, matrix_element):
        """Return the color matrix definition lines for this matrix element. Split rows in chunks of size n."""
        import madgraph.core.color_algebra as color
        if not matrix_element.get('color_matrix'):
            return "\n".join(["      constexpr fptype denom[1] = {1.};", "static const fptype cf[1][1] = {1.};"])
        else:
            color_denominators = matrix_element.get('color_matrix').get_line_denominators()
            denom_string = "      static constexpr fptype denom[ncolor] = {%s};" % ", ".join(["%i" % denom for denom in color_denominators])
            matrix_strings = []
            my_cs = color.ColorString()
            for index, denominator in enumerate(color_denominators):
                # Then write the numerators for the matrix elements
                num_list = matrix_element.get('color_matrix').get_line_numerators(index, denominator)
                matrix_strings.append("{%s}" % ", ".join(["%d" % i for i in num_list]))
            matrix_string = "      static constexpr fptype cf[ncolor][ncolor] = "
            if len( matrix_strings ) > 1 : matrix_string += '{\n      ' + ',\n      '.join(matrix_strings) + '};'
            else: matrix_string += '{' + matrix_strings[0] + '};'
            return "\n".join([denom_string, matrix_string])

    # Replace the export_cpp.OneProcessExporterGPU method (improve formatting)
    def get_initProc_lines(self, matrix_element, color_amplitudes):
        """Get initProc_lines for function definition for gCPPProcess::initProc"""
        initProc_lines = []
        initProc_lines.append("// Set external particle masses for this matrix element")
        for part in matrix_element.get_external_wavefunctions():
            initProc_lines.append("          m_masses.push_back( m_pars->%s );" % part.get('mass'))
        #for i, colamp in enumerate(color_amplitudes):
        #    initProc_lines.append("jamp2_sv[%d] = new double[%d];" % (i, len(colamp))) # N.B. this was commented out already
        return "\n".join(initProc_lines)

    # Replace the export_cpp.OneProcessExporterCPP method (improve formatting)
    def get_helicity_matrix(self, matrix_element):
        """Return the Helicity matrix definition lines for this matrix element"""
        helicity_line = "    , helicities {\n      "; # N.B. SYCL needs access to helicities outside CPPProcess
        helicity_line_list = []
        for helicities in matrix_element.get_helicity_matrix(allow_reverse=True):
            helicity_line_list.append( "{" + ", ".join(['%d'] * len(helicities)) % tuple(helicities) + "}" )
        return helicity_line + ",\n      ".join(helicity_line_list) + "};"

    # Overload the export_cpp.OneProcessExporterGPU method (just to add some comments...)
    def get_reset_jamp_lines(self, color_amplitudes):
        """Get lines to reset jamps"""
        ret_lines = super().get_reset_jamp_lines(color_amplitudes)
        if ret_lines != "" : ret_lines = '    // Reset jamp (reset color flows)\n' + ret_lines # N.B. THIS SHOULD NEVER HAPPEN!
        return ret_lines

#------------------------------------------------------------------------------------

import madgraph.core.helas_objects as helas_objects
import madgraph.iolibs.helas_call_writers as helas_call_writers

# Define a custom HelasCallWriter
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

    # Replace helas_call_writers.GPUFOHelasCallWriter method (improve formatting of gCPPProcess.cu)
    # [GPUFOHelasCallWriter.format_coupling is called by GPUFOHelasCallWriter.get_external_line/generate_helas_call]
    # [GPUFOHelasCallWriter.get_external_line is called by GPUFOHelasCallWriter.get_external]
    # [GPUFOHelasCallWriter.get_external (adding #ifdef CUDA) is called by GPUFOHelasCallWriter.generate_helas_call]
    # [GPUFOHelasCallWriter.generate_helas_call is called by UFOHelasCallWriter.get_wavefunction_call/get_amplitude_call]
    #findcoupling = re.compile('pars->([-]*[\d\w_]+)\s*,')
    def format_coupling(self, call):
        """Format the coupling so any minus signs are put in front"""
        import re
        model = self.get('model')
        if ((not hasattr(self, 'couplings2order')) or (not hasattr(self, 'params2order'))):
            self.couplings2order = {}
            self.params2order = {}
        for coup in re.findall(self.findcoupling, call):
            if coup == 'ZERO':
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
                name = "cIPD"
            else:
                alias = self.couplings2order
                name = "cIPC"
            if coup not in alias:
                alias[coup] = len(alias)
            if name == "cIPD":
                call = call.replace('m_pars->%s%s' % (sign, coup),
                                    '%s%s[%s]' % (sign, name, alias[coup]))
            else:
                call = call.replace('m_pars->%s%s' % (sign, coup),
                                    '%sCOUPs[%s]' % (sign, alias[coup]))
        return call

    # New method for formatting wavefunction/amplitude calls
    # [It would be too complex to modify them in helas_objects.HelasWavefunction/Amplitude.get_call_key]
    @staticmethod
    def format_call(call):
        return call.replace('(','( ').replace(')',' )').replace(',',', ')

    # Replace helas_call_writers.GPUFOHelasCallWriter method (improve formatting)
    def super_get_matrix_element_calls(self, matrix_element, color_amplitudes, multi_channel_map=False):
        """Return a list of strings, corresponding to the Helas calls for the matrix element"""
        import madgraph.core.helas_objects as helas_objects
        import madgraph.loop.loop_helas_objects as loop_helas_objects
        misc.sprint('assert isinstance(matrix_element, helas_objects.HelasMatrixElement)')
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
        res.append('// Reset color flows (reset jamp_sv) at the beginning of a new event or event page')
        res.append('for (size_t i = 0; i < ncolor; i++){ jamp_sv[i] = CXZERO_SV; }')
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
            res.append('\n      // *** DIAGRAM %d OF %d ***' % (diagram.get('number'), len(matrix_element.get('diagrams'))) )
            res.append('\n      // Wavefunction(s) for diagram number %d' % diagram.get('number'))
            res.extend([ self.get_wavefunction_call(wf) for wf in diagram.get('wavefunctions') ]) # N.B. avoid format_call
            if len(diagram.get('wavefunctions')) == 0 : res.append('// (none)')
            if res[-1][-1] == '\n' : res[-1] = res[-1][:-1]
            res.append("\n      // Amplitude(s) for diagram number %d" % diagram.get('number'))
            for amplitude in diagram.get('amplitudes'):
                id_amp +=1
                namp = amplitude.get('number')
                amplitude.set('number', 1)
                res.append(self.get_amplitude_call(amplitude)) # N.B. avoid format_call
                if multi_channel_map: # different code bases #473 (assume this is the same as self.include_multi_channel...)
                    if id_amp in diag_to_config:
                        res.append("#ifdef MGONGPU_SUPPORTS_MULTICHANNEL")
                        res.append("if( channelId == %i ) allNumerators[0] += CXABS2( amp_sv[0] );" % diagram.get('number'))
                        res.append("if( channelId != 0 ) allDenominators[0] += CXABS2( amp_sv[0] );")
                        res.append("#endif")
                else:
                    res.append("#ifdef MGONGPU_SUPPORTS_MULTICHANNEL")
                    res.append("// Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)")
                    res.append("#endif")
                for njamp, coeff in color[namp].items():
                    scoeff = PLUGIN_OneProcessExporter.coeff(*coeff)
                    if scoeff[0] == '+' : scoeff = scoeff[1:]
                    scoeff = scoeff.replace('(','( ')
                    scoeff = scoeff.replace(')',' )')
                    scoeff = scoeff.replace(',',', ')
                    scoeff = scoeff.replace('*',' * ')
                    scoeff = scoeff.replace('/',' / ')
                    if scoeff.startswith('-'): res.append("jamp_sv[%s] -= %samp_sv[0];" % (njamp, scoeff[1:]))
                    else: res.append("jamp_sv[%s] += %samp_sv[0];" % (njamp, scoeff))
            if len(diagram.get('amplitudes')) == 0 : res.append('// (none)')
        #res.append('\n    // *** END OF DIAGRAMS ***' ) # N.B. no longer needed ('COLOR ALGEBRA BELOW')
        return res

    # Overload helas_call_writers.GPUFOHelasCallWriter method (improve formatting)
    def get_matrix_element_calls(self, matrix_element, color_amplitudes, multi_channel_map=False):
        """Return a list of strings, corresponding to the Helas calls for the matrix element"""
        res = self.super_get_matrix_element_calls(matrix_element, color_amplitudes, multi_channel_map)
        for i, item in enumerate(res):
            if item.startswith('# Amplitude'): item='//'+item[1:]
            if not item.startswith('\n') and not item.startswith('#'): res[i]='      '+item
        return res

    # Replace helas_call_writers.GPUFOHelasCallWriter method (improve formatting)
    # [GPUFOHelasCallWriter.format_coupling is called by GPUFOHelasCallWriter.get_external_line/generate_helas_call]
    # [GPUFOHelasCallWriter.get_external_line is called by GPUFOHelasCallWriter.get_external]
    # [=> GPUFOHelasCallWriter.get_external is called by GPUFOHelasCallWriter.generate_helas_call]
    # [GPUFOHelasCallWriter.generate_helas_call is called by UFOHelasCallWriter.get_wavefunction_call/get_amplitude_call]
    first_get_external = True
    def get_external(self, wf, argument):
        misc.sprint('Entering PLUGIN_GPUFOHelasCallWriter.get_external')
        line = self.get_external_line(wf, argument)
        split_line = line.split(',')
        split_line = [ str.lstrip(' ').rstrip(' ') for str in split_line]
        line = ', '.join(split_line)
        if self.first_get_external and ( ( 'mzxxx' in line ) or ( 'pzxxx' in line ) or ( 'xzxxx' in line ) ) :
            self.first_get_external = False
            line2 = line.replace('mzxxx','xxxxx').replace('pzxxx','xxxxx').replace('xzxxx','xxxxx')
            line2 = line2[:line2.find('// NB')]
            split_line2 = line2.split(',')
            split_line2 = [ str.lstrip(' ').rstrip(' ') for str in split_line2]
            split_line2.insert(1, '0') # add parameter fmass=0
            line2 = ', '.join(split_line2)
            text = '#if not defined MGONGPU_TEST_DIVERGENCE\n      %s\n#else\n      if ( ievt %% 2 == 0 )\n        %s\n      else\n        %s\n#endif\n'
            return text % (line, line, line2)
        text = '%s\n'
        return text % line

    # Replace helas_call_writers.GPUFOHelasCallWriter method (vectorize w_sv)
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
            call = call + "( allmomenta + %d,"
            if argument.get('spin') != 1:
                # For non-scalars, need mass and helicity
                call = call + "m_pars->%s, cHel[%d],"
            else:
                # N.B. This seems to be for scalars (spin==1???), pass neither mass nor helicity (#351)
                #call = call + "m_pars->%s,"
                call = call
            call = call + "%+d, w_sv[%d]);"
            if argument.get('spin') == 1:
                misc.sprint("spin")
                misc.sprint(call)
                misc.sprint(
                                (wf.get('number_external')-1,
                                 #wf.get('mass'),
                                 #wf.get('number_external')-1,
                                 # For boson, need initial/final here
                                 (-1) ** (wf.get('state') == 'initial'),
                                 wf.get('me_id')-1))
                return call % \
                                (wf.get('number_external')-1,
                                 #wf.get('mass'),
                                 #wf.get('number_external')-1,
                                 # For boson, need initial/final here
                                 (-1) ** (wf.get('state') == 'initial'),
                                 wf.get('me_id')-1)
            elif argument.is_boson():
                misc.sprint("boson")
                misc.sprint(call)
                misc.sprint( 
                                (wf.get('number_external')-1,
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
                misc.sprint("fermion")
                misc.sprint(call)
                misc.sprint( 
                                (wf.get('number_external')-1,
                                 wf.get('mass'),
                                 wf.get('number_external')-1,
                                 # For fermions, need particle/antiparticle
                                 - (-1) ** wf.get_with_flow('is_part'),
                                 wf.get('me_id')-1))
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
            if wf.get('number_external') == 1 or wf.get('number_external') == 2:
                comment = ' // NB: ' + call + ' only uses pz' # N.B. skip '(not E,px,py)' to avoid interference with comma parsing in get_external
            # Specify namespace for Helas calls
            call = call + '( allmomenta + %d, cHel[%d], %+d, w_sv[%d] );' + comment # vectorize and add comment
            return self.format_coupling(call % \
                                (wf.get('number_external')-1,
                                 wf.get('number_external')-1,
                                 # For fermions, need particle/antiparticle
                                 - (-1) ** wf.get_with_flow('is_part'),
                                 wf.get('me_id')-1))

    # Replace helas_call_writers.GPUFOHelasCallWriter method (vectorize w_sv and amp_sv)
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
        misc.sprint('Entering PLUGIN_GPUFOHelasCallWriter.generate_helas_call')
        if not isinstance(argument, helas_objects.HelasWavefunction) and \
           not isinstance(argument, helas_objects.HelasAmplitude):
            raise self.PhysicsObjectError("get_helas_call must be called with wavefunction or amplitude")
        call = ""
        call_function = None
        if isinstance(argument, helas_objects.HelasAmplitude) and \
           argument.get('interaction_id') == 0:
            misc.sprint('generate_helas_call add amplitude')
            call = "#"
            call_function = lambda amp: call
            self.add_amplitude(argument.get_call_key(), call_function)
            misc.sprint('Leaving PLUGIN_GPUFOHelasCallWriter.generate_helas_call')
            return
        if isinstance(argument, helas_objects.HelasWavefunction) and \
               not argument.get('mothers'):
            misc.sprint('generate_helas_call add wave functions (not mothers)')
            # String is just ixxxxx, oxxxxx, vxxxxx or sxxxxx
            call_function = lambda wf: self.get_external(wf, argument)
        else:
            misc.sprint('generate_helas_call add wave functions')
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
            # (N.B. in the default code these two branches were identical, use a single branch)
            #if isinstance(argument, helas_objects.HelasWavefunction): # e.g. FFV1P0_3 (output is wavefunction)
            #    call = '%(routine_name)s(%(wf)s%(coup)s%(mass)s%(out)s);'
            #else: # e.g. FFV1_0 (output is amplitude)
            #    call = '%(routine_name)s(%(wf)s%(coup)s%(mass)s%(out)s);'
            call = '%(routine_name)s( %(wf)s%(coup)s%(mass)s%(out)s );'
            # compute wf
            arg = {'routine_name': aloha_writers.combine_name('%s' % l[0], l[1:], outgoing, flag, True),
                   'wf': ("w_sv[%%(%d)d], " * len(argument.get('mothers'))) % tuple(range(len(argument.get('mothers')))),
                   'coup': ("m_pars->%%(coup%d)s, " * len(argument.get('coupling'))) % tuple(range(len(argument.get('coupling'))))
                   }
            # determine if this call needs aS-dependent or aS-independent parameters
            usesdepcoupl = None
            for coup in argument.get('coupling'):
                if coup.startswith('-'): coup = coup[1:]
                # Use the same implementation as in UFOModelConverterCPP.prepare_couplings (assume self.model is the same)
                for key, coup_list in self.get('model')['couplings'].items():
                    if coup in coup_list:
                        if "aS" in key:
                            if usesdepcoupl is None: usesdepcoupl = True
                            elif not usesdepcoupl: raise Exception('PANIC! this call seems to use both aS-dependent and aS-independent couplings?')
                        else:
                            if usesdepcoupl is None: usesdepcoupl = False
                            elif usesdepcoupl: raise Exception('PANIC! this call seems to use both aS-dependent and aS-independent couplings?')
            if usesdepcoupl is None: raise Exception('PANIC! could not determine if this call uses aS-dependent or aS-independent couplings?')
            if isinstance(argument, helas_objects.HelasWavefunction):
                arg['out'] = 'w_sv[%(out)d]'
                if aloha.complex_mass:
                    arg['mass'] = "m_pars->%(CM)s, "
                else:
                    arg['mass'] = "m_pars->%(M)s, m_pars->%(W)s, "
            else:
                arg['out'] = '&amp_sv[%(out)d]'
                arg['out2'] = 'amp_sv[%(out)d]'
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
        misc.sprint('Leaving PLUGIN_GPUFOHelasCallWriter.generate_helas_call')

# Use the custom HelasCallWriter
DEFAULT_GPUFOHelasCallWriter = helas_call_writers.GPUFOHelasCallWriter
helas_call_writers.GPUFOHelasCallWriter = PLUGIN_GPUFOHelasCallWriter

#------------------------------------------------------------------------------------
