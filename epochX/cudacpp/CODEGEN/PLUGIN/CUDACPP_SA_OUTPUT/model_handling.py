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

# AV - replace aloha_writers.Declaration_list.is_used (disable caching to be on the safe side)
# (NB class Declaration_list(set) is a set of (type, name) pairs!)
def PLUGIN_Declaration_list_is_used(self, var):
    ###if hasattr(self, 'var_name'): return var in self.var_name # AV why was this needed? disable caching to be on the safe side
    self.var_name = [name for type,name in self]
    return var in self.var_name

DEFAULT_Declaration_list_is_used = aloha_writers.Declaration_list.is_used
aloha_writers.Declaration_list.is_used = PLUGIN_Declaration_list_is_used

#------------------------------------------------------------------------------------

# AV - decorate aloha_writers.Declaration_list.add (add optional debug printout)
def PLUGIN_Declaration_list_add(self, obj):
    #print( 'ADDING ', obj) # FOR DEBUGGING
    #assert( obj[1] != 'P3' ) # FOR DEBUGGING (check MG5_debug to see where OM3, TMP3, P3 etc were added)
    return DEFAULT_Declaration_list_add(self, obj)

DEFAULT_Declaration_list_add = aloha_writers.Declaration_list.add
aloha_writers.Declaration_list.add = PLUGIN_Declaration_list_add

#------------------------------------------------------------------------------------

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

    # AV - add vector types
    type2def['double_v'] = 'fptype_sv'
    type2def['complex_v'] = 'cxtype_sv'

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
            elif number == -1: out = '- one' # AV
            elif number == 2: out = 'two' # AV
            elif number == -2: out = '- two' # AV
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
                        out = '- cI'
                    else: 
                        out = '%s * cI' % self.change_number_format(number.imag)
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
                type = self.type2def[format[5:]] # double or complex (instead of list_double or list_complex)
                if not argname.startswith('COUP'): type = self.type2def[format[5:]+'_v'] # AV vectorize (double_v or complex_v)
                list_arg = '[]'
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
            output = '%(doublec)s%(pointer_vertex)s vertex' % { # AV vectorize
                'doublec':self.type2def['complex_v'],
                'pointer_vertex': self.type2def['pointer_vertex']}
            comment_output = 'amplitude \'vertex\''
        else:
            ###output = '%(doublec)s %(spin)s%(id)d[]' % {
            output = '%(doublec)s %(spin)s%(id)d[]' % { # AV vectorize
                     'doublec': self.type2def['complex_v'],
                     'spin': self.particles[self.outgoing -1],
                     'id': self.outgoing}
            ###self.declaration.add(('list_complex', output)) # AV BUG FIX - THIS IS NOT NEEDED AND IS WRONG (adds name 'cxtype_sv V3[]')
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
        out.write('    mgDebug( 0, __FUNCTION__ );\n') # AV
        argument_var = [name for type,name in self.call_arg]
        # define the complex number CI = 0+1j
        if add_i:
            ###out.write(self.ci_definition)
            out.write('    ' + self.ci_definition) # AV
        codedict = {} # AV allow delayed declaration with initialisation
        for type, name in self.declaration.tolist():
            ###print(name) # FOR DEBUGGING
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
                codedict[fullname] = '%s %s' % (self.type2def[type+"_v"], fullname) # AV vectorize, add to codedict
            ###print(fullname, codedict[fullname]) # FOR DEBUGGING
            if self.nodeclare:
                self.declaration.codedict = codedict # AV new behaviour (delayed declaration with initialisation)
            else:
                out.write('    %s;\n' % codedict[fullname] ) # AV old behaviour (separate declaration with no initialization)
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
                ###out.write("    OM{0} = ( M{0} != {1} ? {2} / ( M{0} * M{0} ) : {1} );\n".format( # AV use ternary in OM3
                ###    ###i+1, self.change_number_format(0), self.change_number_format(1)))
                ###    i+1, '0.', '1.')) # AV force scalar "1." instead of vector "one"
                declname = 'OM%s' % (i+1) # AV
                if self.nodeclare: declname = 'const ' + self.declaration.codedict[declname] # AV
                out.write("    {3} = ( M{0} != {1} ? {2} / ( M{0} * M{0} ) : {1} );\n".format( # AV use ternary in OM3
                    i+1, '0.', '1.', declname)) # AV force scalar "1." instead of vector "one", add declaration
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

    # AV - modify aloha_writers.ALOHAWriterForCPP method (improve formatting, add delayed declaration with initialisation)
    # This affects 'P1[0] = ' in HelAmps_sm.cu
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
                if self.nodeclare:
                    mydict['declnamedenom'] = 'const %s denom' % self.type2def['complex_v'] # AV
                else:
                    mydict['declnamedenom'] = 'denom' # AV
                    self.declaration.add(('complex','denom'))
                if not aloha.complex_mass:
                    # This affects 'denom = COUP' in HelAmps_sm.cu
                    if self.routine.denominator:
                        out.write('    %(declnamedenom)s = %(pre_coup)s%(coup)s%(post_coup)s / (%(denom)s)\n' % mydict) # AV
                    else:
                        out.write('    %(declnamedenom)s = %(pre_coup)s%(coup)s%(post_coup)s / ( (P%(i)s[0] * P%(i)s[0] ) - ( P%(i)s[1] * P%(i)s[1] ) - ( P%(i)s[2] * P%(i)s[2] ) - ( P%(i)s[3] * P%(i)s[3] ) - M%(i)s * ( M%(i)s - cI * W%(i)s ) );\n' % mydict) # AV
                else:
                    if self.routine.denominator:
                        raise Exception('modify denominator are not compatible with complex mass scheme')                
                    # This affects 'denom = COUP' in HelAmps_sm.cu
                    out.write('    %(declnamedenom)s = %(pre_coup)s%(coup)s%(post_coup)s / ( (P%(i)s[0] * P%(i)s[0] ) - ( P%(i)s[1] *P%(i)s[1] ) - ( P%(i)s[2] * P%(i)s[2] ) - ( P%(i)s[3] * P%(i)s[3] ) - ( M%(i)s * M%(i)s ) );\n' % mydict) # AV
                ###self.declaration.add(('complex','denom')) # AV moved earlier (or simply removed)
                if aloha.loop_mode: ptype = 'list_complex'
                else: ptype = 'list_double'
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

    # AV - modify aloha_writers.WriteALOHA method (improve formatting)
    # This affects 'V1[2] = ' and 'F1[2] = ' in HelAmps_sm.cu
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
            ###file_str.write('*(')
            file_str.write(' * ( ') # AV
        else:
            ###file_str.write('(')
            file_str.write('( ') # AV
        ###print('."'+file_str.getvalue()+'"') # AV - FOR DEBUGGING
        first=True
        for value, obj_list in data.items():
            ###print('.."' + str(value) + '" "' + str(obj_list) + '"') # AV - FOR DEBUGGING
            add= ' + '
            if value not in  [-1,1]:
                nb_str = self.change_number_format(value)
                if nb_str[0] in ['+','-']:
                    ###file_str.write(' '+nb_str) # AV
                    file_str.write(nb_str) # AV
                else:
                    ###file_str.write('+')
                    file_str.write('+' if first else ' + ') # AV
                    file_str.write(nb_str)
                ###file_str.write('*(')
                file_str.write(' * ( ') # AV (eg '+ cI * (V3[4])')
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
            # AV comment: write_obj here also adds calls declaration_add (via change_var_format) - example: OM3
            ###print('..."'+file_str.getvalue()+'"') # AV - FOR DEBUGGING
            file_str.write(add.join([self.write_obj(obj, prefactor=False) for obj in obj_list]))
            if value not in [1,-1]:
                ###file_str.write(')')
                file_str.write(' )') # AV
        if number:
            total = sum(number)
            file_str.write('+ %s' % self.change_number_format(total))
        ###file_str.write(')')
        file_str.write(' )') # AV
        ###print('...."'+file_str.getvalue()+'"') # AV - FOR DEBUGGING
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
        if len(res_strings) == 0 : res_strings.append('// (none)')
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
        replace_dict['hel_amps_cc'] = "#include \"HelAmps_%s.cc\"" % self.model_name # AV
        coupling = [''] * len(self.couplings2order)
        params = [''] * len(self.params2order)
        for coup, pos in self.couplings2order.items():
            coupling[pos] = coup
        ###coup_str = "static cxtype tIPC[%s] = {pars->%s};\n"\
        ###    %(len(self.couplings2order), ',pars->'.join(coupling))
        coup_str = "const cxtype tIPC[%s] = { cxmake( m_pars->%s ) };\n"\
            %(len(self.couplings2order), ' ), cxmake( m_pars->'.join(coupling)) # AV
        for para, pos in self.params2order.items():
            params[pos] = para
        ###param_str = "static double tIPD[%s] = {pars->%s};\n"\
        ###    %(len(self.params2order), ',pars->'.join(params))
        param_str = "    const fptype tIPD[%s] = { (fptype)m_pars->%s };"\
            %(len(self.params2order), ', (fptype)m_pars->'.join(params)) # AV
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
            ret_lines.append("""
  // Evaluate |M|^2 for each subprocess
  // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event(s)
  __device__
  INLINE
  void calculate_wavefunctions( int ihel,
                                const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                                fptype_sv* allMEs            // output: allMEs[npagM][neppM], final |M|^2 averaged over helicities
#ifndef __CUDACC__
                                , const int nevt             // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
                                )
  //ALWAYS_INLINE // attributes are not permitted in a function definition
  {
    using namespace MG5_sm;
    mgDebug( 0, __FUNCTION__ );
#ifndef __CUDACC__
    //printf( "calculate_wavefunctions: nevt %d\\n", nevt );
#endif\n""")
            ret_lines.append("    // The number of colors")
            ret_lines.append("    constexpr int ncolor = %i;" % len(color_amplitudes[0]))
            ret_lines.append("""
    // Local TEMPORARY variables for a subset of Feynman diagrams in the given CUDA event (ievt) or C++ event page (ipagV)
    // [NB these variables are reused several times (and re-initialised each time) within the same event or event page]
    cxtype_sv w_sv[nwf][nw6]; // particle wavefunctions within Feynman diagrams (nw6 is often 6, the dimension of spin 1/2 or spin 1 particles)
    cxtype_sv amp_sv[1]; // invariant amplitude for one given Feynman diagram

    // Local variables for the given CUDA event (ievt) or C++ event page (ipagV)
    cxtype_sv jamp_sv[ncolor] = {}; // sum of the invariant amplitudes for all Feynman diagrams in the event or event page

    // === Calculate wavefunctions and amplitudes for all diagrams in all processes - Loop over nevt events ===
#ifndef __CUDACC__
    const int npagV = nevt / neppV;
    // ** START LOOP ON IPAGV **
#ifdef _OPENMP
    // (NB gcc9 or higher, or clang, is required)
    // - default(none): no variables are shared by default
    // - shared: as the name says
    // - private: give each thread its own copy, without initialising
    // - firstprivate: give each thread its own copy, and initialise with value from outside
#pragma omp parallel for default(none) shared(allmomenta,allMEs,cHel,cIPC,cIPD,ihel,npagV) private (amp_sv,w_sv,jamp_sv)
#endif
    for ( int ipagV = 0; ipagV < npagV; ++ipagV )
#endif""")
            ret_lines.append('    {') # NB This is closed in process_matrix.inc
            helas_calls = self.helas_call_writer.get_matrix_element_calls(\
                                                    self.matrix_elements[0],
                                                    color_amplitudes[0]
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
        self.edit_processidfile() # AV new file (NB this is Sigma-specific, should not be a symlink to Subprocesses)
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
        # AV May remove replace_dict as no replacement is done in check_sa.cc (in upstream Madgraph)
        ###replace_dict = {}
        ###replace_dict['nexternal'], _ = self.matrix_elements[0].get_nexternal_ninitial()
        ###replace_dict['model'] = self.model_name
        ###replace_dict['numproc'] = len(self.matrix_elements)
        ff = open(pjoin(self.path, 'check_sa.cc'),'w')
        ff.write(template)
        ###ff.write(template % replace_dict) # AV normally this should be used! (and % should be %% in check_sa.cc)
        ff.close()

    # AV - add debug printouts over the export_cpp.OneProcessExporterGPU method
    def edit_mgonGPU(self):
        """Generate mgOnGpuConfig.h"""
        misc.sprint('Entering PLUGIN_OneProcessExporter.edit_mgonGPU')
        ###misc.sprint('  template_path=%s'%self.template_path) # look for gpu/mgOnGpuConfig.h here
        return super().edit_mgonGPU()

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
        replace_dict['hel_amps_h'] = "#include \"HelAmps_%s.h\"" % self.model_name # AV
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
                return '-' # AV keep default (eg jamp_sv[0] += -amp_sv[0])
        assert(False)
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
            ###return "\n".join(["static const double denom[1] = {1.};", "static const double cf[1][1] = {1.};"])
            return "\n".join(["      static constexpr fptype denom[1] = {1.};", "static const fptype cf[1][1] = {1.};"]) # AV
        else:
            color_denominators = matrix_element.get('color_matrix').\
                                                 get_line_denominators()
            ###denom_string = "static const double denom[ncolor] = {%s};" % ",".join(["%i" % denom for denom in color_denominators])
            denom_string = "      static constexpr fptype denom[ncolor] = {%s};" % ", ".join(["%i" % denom for denom in color_denominators]) # AV
            matrix_strings = []
            my_cs = color.ColorString()
            for index, denominator in enumerate(color_denominators):
                # Then write the numerators for the matrix elements
                num_list = matrix_element.get('color_matrix').get_line_numerators(index, denominator)
                ###matrix_strings.append("{%s}" % ",".join(["%d" % i for i in num_list]))
                matrix_strings.append("{%s}" % ", ".join(["%d" % i for i in num_list])) # AV
            ###matrix_string = "static const double cf[ncolor][ncolor] = {" + ",".join(matrix_strings) + "};"
            matrix_string = "      static constexpr fptype cf[ncolor][ncolor] = " # AV
            if len( matrix_strings ) > 1 : matrix_string += '{\n      ' + ',\n      '.join(matrix_strings) + '};' # AV
            else: matrix_string += '{' + matrix_strings[0] + '};' # AV
            return "\n".join([denom_string, matrix_string])

    # AV - replace the export_cpp.OneProcessExporterGPU method (improve formatting)
    def get_initProc_lines(self, matrix_element, color_amplitudes):
        """Get initProc_lines for function definition for gCPPProcess::initProc"""
        initProc_lines = []
        initProc_lines.append("// Set external particle masses for this matrix element")
        for part in matrix_element.get_external_wavefunctions():
            ###initProc_lines.append("mME.push_back(pars->%s);" % part.get('mass'))
            initProc_lines.append("    m_masses.push_back( m_pars->%s );" % part.get('mass')) # AV
        ###for i, colamp in enumerate(color_amplitudes):
        ###    initProc_lines.append("jamp2_sv[%d] = new double[%d];" % (i, len(colamp))) # AV - this was commented out already
        return "\n".join(initProc_lines)

    # AV - replace the export_cpp.OneProcessExporterCPP method (improve formatting)
    def get_helicity_matrix(self, matrix_element):
        """Return the Helicity matrix definition lines for this matrix element"""
        ###helicity_line = "static const int helicities[ncomb][nexternal] = {";
        helicity_line = "    static constexpr short helicities[ncomb][mgOnGpu::npar] = {\n      "; # AV (this is tHel)
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

import madgraph.core.helas_objects as helas_objects
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
                name = "cIPD"
            else: 
                alias = self.couplings2order
                name = "cIPC"
            if coup not in alias:
                alias[coup] = len(alias)
            if name == "cIPD":
                ###call = call.replace('pars->%s%s' % (sign, coup), 
                call = call.replace('m_pars->%s%s' % (sign, coup), # AV
                                    '%s%s[%s]' % (sign, name, alias[coup]))
            else:
                ###call = call.replace('pars->%s%s' % (sign, coup), 
                call = call.replace('m_pars->%s%s' % (sign, coup), # AV
                                    ###'%scxtype(cIPC[%s],cIPC[%s])' % 
                                    '%scxmake( cIPC[%s], cIPC[%s] )' % 
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
        res.append('// Reset color flows (reset jamp_sv) at the beginning of a new event or event page')
        res.append('for( int i=0; i<ncolor; i++ ){ jamp_sv[i] = cxzero_sv(); }')
        for diagram in matrix_element.get('diagrams'):
            ###print('DIAGRAM %3d: #wavefunctions=%3d, #diagrams=%3d' %
            ###      (diagram.get('number'), len(diagram.get('wavefunctions')), len(diagram.get('amplitudes')) )) # AV - FOR DEBUGGING
            res.append('\n      // *** DIAGRAM %d OF %d ***' % (diagram.get('number'), len(matrix_element.get('diagrams'))) ) # AV
            res.append('\n      // Wavefunction(s) for diagram number %d' % diagram.get('number')) # AV
            res.extend([ self.get_wavefunction_call(wf) for wf in diagram.get('wavefunctions') ]) # AV new: avoid format_call
            ###res.extend([ self.format_call(self.get_wavefunction_call(wf)) for wf in diagram.get('wavefunctions') ]) # AV old: use format_call
            if len(diagram.get('wavefunctions')) == 0 : res.append('// (none)') # AV
            ###res.append("# Amplitude(s) for diagram number %d" % diagram.get('number'))
            res.append("\n      // Amplitude(s) for diagram number %d" % diagram.get('number'))
            for amplitude in diagram.get('amplitudes'):
                namp = amplitude.get('number')
                amplitude.set('number', 1)
                res.append(self.get_amplitude_call(amplitude)) # AV new: avoid format_call
                ###res.append(self.format_call(self.get_amplitude_call(amplitude))) # AV old: use format_call
                for njamp, coeff in color[namp].items():
                    ###res.append("jamp[%s] += %samp[0];" % (njamp, export_cpp.OneProcessExporterGPU.coeff(*coeff)))
                    ###res.append("jamp[%s] += %samp[0];" % (njamp, PLUGIN_OneProcessExporter.coeff(*coeff)))
                    ###res.append("jamp_sv[%s] += %samp_sv[0];" % (njamp, PLUGIN_OneProcessExporter.coeff(*coeff))) # AV vectorize
                    scoeff = PLUGIN_OneProcessExporter.coeff(*coeff) # AV
                    if scoeff.startswith('-'): res.append("jamp_sv[%s] -= %samp_sv[0];" % (njamp, scoeff[1:])) # AV
                    else: res.append("jamp_sv[%s] += %samp_sv[0];" % (njamp, scoeff)) # AV
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
            if not item.startswith('\n') and not item.startswith('#'): res[i]='      '+item
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
        line = ', '.join(split_line) # AV (for CUDA)
        ###split_line.insert(-1, ' ievt')
        split_line.insert(-1, 'ipagV') # AV (for C++)
        if self.first_get_external and ( ( 'mzxxx' in line ) or ( 'pzxxx' in line ) or ( 'xzxxx' in line ) ) :
            self.first_get_external = False
            line2 = line.replace('mzxxx','xxxxx').replace('pzxxx','xxxxx').replace('xzxxx','xxxxx')
            line2 = line2[:line2.find('// NB')]
            split_line2 = line2.split(',')
            split_line2 = [ str.lstrip(' ').rstrip(' ') for str in split_line2] # AV
            split_line2.insert(1, '0') # add parameter fmass=0
            line2 = ', '.join(split_line2)
            text = '#ifdef __CUDACC__\n#ifndef MGONGPU_TEST_DIVERGENCE\n      %s\n#else\n      if ( ( blockDim.x * blockIdx.x + threadIdx.x ) %% 2 == 0 )\n        %s\n      else\n        %s\n#endif\n#else\n      %s\n#endif\n' # AV
            return text % (line, line, line2, ', '.join(split_line))
        ###text = '\n#ifdef __CUDACC__\n    %s    \n#else\n    %s\n#endif \n'
        text = '#ifdef __CUDACC__\n      %s\n#else\n      %s\n#endif\n' # AV
        return text % (line, ', '.join(split_line))
    
    # AV - replace helas_call_writers.GPUFOHelasCallWriter method (vectorize w_sv)
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
            call = call + "(allmomenta,"
            if argument.get('spin') != 1:
                # For non-scalars, need mass and helicity
                call = call + "m_pars->%s, cHel[ihel][%d],"
            else:
                call = call + "m_pars->%s,"
            ###call = call + "%+d,w[%d], %d);"
            call = call + "%+d, w_sv[%d], %d);" # AV vectorize
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
                comment = '' # AV
            call = call + 'x' * (6 - len(call))
            if wf.get('number_external') == 1 or wf.get('number_external') == 2: # AV
                comment = ' // NB: ' + call + ' only uses pz' # AV skip '(not E,px,py)' to avoid interference with comma parsing in get_external
            # Specify namespace for Helas calls
            ###call = call + "(allmomenta, cHel[ihel][%d],%+d,w[%d],%d);"
            call = call + '( allmomenta, cHel[ihel][%d], %+d, w_sv[%d], %d );' + comment # AV vectorize and add comment
            return self.format_coupling(call % \
                                (wf.get('number_external')-1,
                                 # For fermions, need particle/antiparticle
                                 - (-1) ** wf.get_with_flow('is_part'),
                                 wf.get('me_id')-1,
                                 wf.get('number_external')-1))

    # AV - replace helas_call_writers.GPUFOHelasCallWriter method (vectorize w_sv and amp_sv)
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
            raise self.PhysicsObjectError("get_helas_call must be called with wavefunction or amplitude")
        call = ""
        call_function = None
        if isinstance(argument, helas_objects.HelasAmplitude) and \
           argument.get('interaction_id') == 0:
            call = "#"
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
            arg = {'routine_name': aloha_writers.combine_name(\
                                            '%s' % l[0], l[1:], outgoing, flag,True),
                   ###'wf': ("w[%%(%d)d]," * len(argument.get('mothers'))) % tuple(range(len(argument.get('mothers')))),
                   'wf': ("w_sv[%%(%d)d], " * len(argument.get('mothers'))) % tuple(range(len(argument.get('mothers')))), # AV
                   ###'coup': ("pars->%%(coup%d)s," * len(argument.get('coupling'))) % tuple(range(len(argument.get('coupling'))))
                   'coup': ("m_pars->%%(coup%d)s, " * len(argument.get('coupling'))) % tuple(range(len(argument.get('coupling')))) # AV
                   }
            if isinstance(argument, helas_objects.HelasWavefunction):
                ###arg['out'] = 'w[%(out)d]'
                arg['out'] = 'w_sv[%(out)d]'
                if aloha.complex_mass:
                    ###arg['mass'] = "pars->%(CM)s,"
                    arg['mass'] = "m_pars->%(CM)s, " # AV
                else:
                    ###arg['mass'] = "pars->%(M)s,pars->%(W)s,"
                    arg['mass'] = "m_pars->%(M)s, m_pars->%(W)s, " # AV
            else:        
                ###arg['out'] = '&amp[%(out)d]'
                ###arg['out2'] = 'amp[%(out)d]'
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

# AV - use the custom HelasCallWriter
DEFAULT_GPUFOHelasCallWriter = helas_call_writers.GPUFOHelasCallWriter
helas_call_writers.GPUFOHelasCallWriter = PLUGIN_GPUFOHelasCallWriter

#------------------------------------------------------------------------------------
