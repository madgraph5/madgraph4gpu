# Copyright (C) 2020-2025 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Oct 2026) for the MadVFold project (arXiv:2609.23500).
#
# AV vfold (step 9i+9j): code generation of the per-event "Born extras" block
# [colour-linked Borns, amp2, jamp2] that the cudacpp library returns to the MadFKS
# Fortran SBORN (fbridgegetbsf_ and fbridgegetbornextras_), so that mode C of vfold
# no longer needs to run the scalar Fortran SBORN at all.
#
# The colour-linked Borns B_mn = <M|T_m.T_n|M> are those of the MadFKS Fortran routines
# B_SF_001..B_SF_NNN, computed here with the same MG5 colour algebra as in export_fks.py
# (fks_common.find_color_links and insert_color_links). Each B_SF computes
# Re(sum_ij CF(j,i) JAMP2(j) conj(JAMP1(i)))/IDEN from the Born colour flows JAMP1 and the
# colour-linked flows JAMP2. As the M = sum_l c_l JAMP1(l) is linear in the Born colour basis
# c_l and JAMP2 is obtained by inserting the link in the colour factors, each link is the
# quadratic form Q_il = <c_i|T_m.T_n|c_l> in the Born jamps, i.e. Q = CF D where D(j,l) is
# the component of T_m.T_n c_l on the j-th element of the linked colour basis. This is
# derived here with exact rational arithmetic (by inserting the link in the Born colour basis
# elements themselves) and is hardcoded in CPPProcess.cc.

from fractions import Fraction

import madgraph.core.color_algebra as color_algebra
import madgraph.fks.fks_common as fks_common

NUMBER_WORDS = ['zero', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine', 'ten',
                'eleven', 'twelve', 'thirteen', 'fourteen', 'fifteen', 'sixteen', 'seventeen', 'eighteen',
                'nineteen', 'twenty']

class VfoldBornExtrasError(Exception):
    pass

#------------------------------------------------------------------------------------

def is_vfold(output_options):
    """Is the --vfold output option set? (e.g. 'output madevent_simd PROC.mad --vfold')"""
    vfold = output_options.get('vfold', False)
    if isinstance(vfold, str): vfold = vfold.lower() in ('true', 't', '1', 'yes', 'on')
    return bool(vfold)

def filter_vfold_blocks(text, vfold):
    """Keep or remove the vfold-only blocks of a template file: the lines between '%%VFOLD_BEGIN'
    and '%%VFOLD_END' are kept only with --vfold, those between an optional '%%VFOLD_ELSE' and
    '%%VFOLD_END' only without --vfold (the marker lines themselves are always removed)"""
    if '%%VFOLD_' not in text: return text
    out = []
    block = None # None (outside a vfold block), 'BEGIN' or 'ELSE'
    for line in text.splitlines(keepends=True):
        marker = line.strip()
        if marker == '%%VFOLD_BEGIN' and block is None: block = 'BEGIN'
        elif marker == '%%VFOLD_ELSE' and block == 'BEGIN': block = 'ELSE'
        elif marker == '%%VFOLD_END' and block is not None: block = None
        elif marker.startswith('%%VFOLD_'): raise VfoldBornExtrasError('unexpected %s in a template file' % marker)
        elif block is None or ( block == 'BEGIN' ) == vfold: out.append(line)
    if block is not None: raise VfoldBornExtrasError('missing %%VFOLD_END in a template file')
    return ''.join(out)

#------------------------------------------------------------------------------------

def number_word(n):
    """Spell out small numbers (e.g. 'ten'), keep digits for larger ones"""
    return NUMBER_WORDS[n] if n < len(NUMBER_WORDS) else '%d' % n

def cxadd(a, b):
    return (a[0] + b[0], a[1] + b[1])

def cxmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])

CXZERO = (Fraction(0), Fraction(0))

def color_coeff(coeff, Nc_value=3):
    """Complex Fraction (re, im) for a colour amplitude coefficient (ff_number, frac, is_imaginary, Nc_power)"""
    ff_number, frac, is_imaginary, Nc_power = coeff
    total = Fraction(ff_number) * Fraction(frac) * Fraction(Nc_value) ** Nc_power
    return (Fraction(0), total) if is_imaginary else (total, Fraction(0))

def jamp_coefficients(color_amplitudes, namp):
    """Matrix F[icol][iamp] such that jamp_icol = sum_iamp F[icol][iamp] * amp_iamp"""
    F = [[CXZERO] * namp for _ in color_amplitudes]
    for icol, col in enumerate(color_amplitudes):
        for coeff, iamp in col:
            F[icol][iamp - 1] = cxadd(F[icol][iamp - 1], color_coeff(coeff))
    return F

def format_fraction(f):
    """C++ literal for an exact fraction (e.g. '6' or '-4. / 3.')"""
    if f.denominator == 1: return '%d' % f.numerator
    return '%d. / %d.' % (f.numerator, f.denominator)

def format_ranges(numbers):
    """Compact list of numbers (e.g. '1-4,13-16')"""
    ranges = []
    for n in numbers:
        if ranges and n == ranges[-1][1] + 1: ranges[-1][1] = n
        else: ranges.append([n, n])
    return ','.join(['%d' % a if a == b else '%d-%d' % (a, b) for a, b in ranges])

def join_words(items):
    """Join a list of strings as 'a', 'a and b' or 'a, b and c'"""
    if len(items) <= 1: return ''.join(items)
    return ', '.join(items[:-1]) + ' and ' + items[-1]

#------------------------------------------------------------------------------------

def get_amp2_lines(iamp):
    """The AMP2 update after the jamp updates of amplitude iamp (starting at 0)
    [NB: step 14i: the Fortran AMP2 has one entry per amplitude, not per diagram (a diagram with
    a 4-gluon vertex has several amplitudes), and iamp is the index of the Fortran AMP]
    [NB: in CUDA the one-helicity buffer is not reset, so its value is assigned]"""
    return ['#ifndef MGONGPUCPP_GPUIMPL',
            '// AV vfold: Fortran AMP2(%d)' % (iamp + 1),
            'if( bsfx_sv ) bsfx_sv[nbsf + %d] += cxabs2( amp_sv[0] );' % (iamp),
            '#else',
            '// AV vfold: Fortran AMP2(%d) for this helicity' % (iamp + 1),
            'if( amp2_ev ) amp2_ev[%d * nevt] = cxabs2( amp_sv[0] );' % (iamp),
            '#endif']

#------------------------------------------------------------------------------------

class VfoldBornExtras(object):
    """The colour-linked Borns, amp2 and jamp2 of one Born matrix element (vfold step 9i+9j)"""

    def __init__(self, matrix_element):
        self.matrix_element = matrix_element
        amplitude = matrix_element.get('base_amplitude')
        process = amplitude.get('process')
        self.process_string = process.base_string()
        self.namp = matrix_element.get_number_of_amplitudes()
        col_basis = matrix_element.get('color_basis')
        color_amplitudes = matrix_element.get_color_amplitudes()
        self.ncolor = len(color_amplitudes)
        self.denominator = matrix_element.get_denominator_factor()
        F = jamp_coefficients(color_amplitudes, self.namp)
        # Diagram groups with identical colour coefficients (for the comments)
        self.diagram_groups = self.get_diagram_groups(F)
        # The colour matrix of the Born (for the comments)
        self.cf = self.get_color_matrix(matrix_element.get('color_matrix'))
        # The colour links, in the same order as the MadFKS B_SF_00N
        legs = process.get('legs')
        model = process.get('model')
        links = fks_common.find_color_links(fks_common.to_fks_legs(legs, model), symm=True, pert='QCD')
        self.colored_legs = sorted(set([leg.get('number') for link in links for leg in link['legs']]))
        # Insert the links in the Born colour basis elements c_l themselves (one "diagram" per element)
        basis_elements = []
        for immutable in sorted(col_basis.keys()): # same order as the jamps and the colour matrix
            col_str = color_algebra.ColorString()
            col_str.from_immutable(immutable)
            basis_elements.append({(0,): col_str})
        color_links = fks_common.insert_color_links(col_basis, basis_elements, links)
        self.links = []
        self.qforms = []
        for link in color_links:
            link_basis = link['link_basis']
            cf2 = link['link_matrix'].col_matrix_fixed_Nc
            # D(j,l): component of T_m.T_n c_l on the j-th element of the linked colour basis
            D = [[CXZERO] * self.ncolor for _ in link_basis]
            for j, immutable in enumerate(sorted(link_basis.keys())):
                for (l, col_chain, coeff, is_imaginary, Nc_power, loop_Nc_power) in link_basis[immutable]:
                    D[j][l] = cxadd(D[j][l], color_coeff((1, coeff, is_imaginary, Nc_power)))
            Q = []
            for i in range(self.ncolor):
                Qrow = []
                for l in range(self.ncolor):
                    s = CXZERO
                    for j in range(len(link_basis)):
                        value, is_imaginary = cf2[(i, j)][0], cf2[(i, j)][1]
                        cfij = (Fraction(0), value) if is_imaginary else (value, Fraction(0))
                        s = cxadd(s, cxmul(cfij, D[j][l]))
                    Qrow.append(s)
                Q.append(Qrow)
            # Re(J^+ Q J) only involves the hermitian part of Q, which must be real (symmetric)
            for i in range(self.ncolor):
                for l in range(self.ncolor):
                    if Q[i][l][1] != -Q[l][i][1]:
                        raise VfoldBornExtrasError('complex colour-linked Born quadratic form for link %s' % link['link'])
            self.links.append(link['link'])
            self.qforms.append([[(Q[i][l][0] + Q[l][i][0]) / 2 for l in range(self.ncolor)] for i in range(self.ncolor)])
        self.nbsf = len(self.links)
        # Step 14g: the second component ANS(2) of the Fortran SBORN, for each possible massless mother
        # of an FKS splitting: the final-state massless coloured legs, in the order of the Born legs
        # (the same legs as in vfold_ans2_index in vfold_bsf.inc, from madgraph/iolibs/export_fks.py).
        # ANS(2) is the interference of the helicities of this leg: for each helicity with the helicity
        # of this leg as in the first combination (BACK_HEL in Fortran), cAns2Partner gives the helicity
        # that only differs by the opposite helicity of this leg
        self.ans2legs = [leg.get('number') for leg in legs if leg.get('state')
                         and model.get_particle(leg.get('id')).get('mass').lower() == 'zero'
                         and model.get_particle(leg.get('id')).get('color') != 1]
        self.nans2 = len(self.ans2legs)
        helicities = [tuple(hel) for hel in matrix_element.get_helicity_matrix(allow_reverse=True)]
        self.ncomb = len(helicities)
        self.ans2partner = []
        for number in self.ans2legs:
            ileg = number - 1
            partner = []
            for hel in helicities:
                if hel[ileg] != helicities[0][ileg]: partner.append(-1)
                else: partner.append(helicities.index(hel[:ileg] + (-hel[ileg],) + hel[ileg + 1:]))
            self.ans2partner.append(partner)
        self.nbx = self.nbsf + self.namp + self.ncolor + 2 * self.nans2
        # The terms of the quadratic forms: diagonal first, then the upper triangle (skip terms that are always zero)
        self.terms = [(i, i) for i in range(self.ncolor)]
        self.terms += [(i, l) for i in range(self.ncolor) for l in range(i + 1, self.ncolor)
                       if any(q[i][l] != 0 for q in self.qforms)]
        self.all_integer = all(q[i][l].denominator == 1 for q in self.qforms for (i, l) in self.terms)

    @staticmethod
    def get_diagram_groups(F):
        groups = {}
        for iamp in range(len(F[0]) if F else 0):
            key = tuple(row[iamp] for row in F)
            groups.setdefault(key, []).append(iamp + 1)
        return sorted(groups.values())

    @staticmethod
    def get_color_matrix(color_matrix):
        n = len(color_matrix._col_basis1) if color_matrix else 0
        return [[color_matrix.col_matrix_fixed_Nc[(i, j)][0] for j in range(n)] for i in range(n)]

    #----------------------------------------------------------------------------------

    def suffix(self, i, l):
        """Suffix of the names of a term (e.g. '12' for jamp0*jamp1)"""
        if self.ncolor < 10: return '%d%d' % (i + 1, l + 1)
        return '%d_%d' % (i + 1, l + 1)

    def term_formula(self, i, l):
        if i == l: return 'S%s[k]*|jamp%d|^2' % (self.suffix(i, l), i)
        return '2*S%s[k]*Re(conj(jamp%d)*jamp%d)' % (self.suffix(i, l), i, l)

    def cf_string(self):
        return '{' + ','.join(['{' + ','.join([str(x) for x in row]) + '}' for row in self.cf]) + '}'

    def coefficient_arrays(self, indent):
        """The static constexpr arrays of the coefficients of all links, one per term"""
        lines = []
        for (i, l) in self.terms:
            values = ', '.join([format_fraction(q[i][l]) for q in self.qforms])
            lines.append('%sstatic constexpr fptype cBsfS%s[nbsf] = { %s };' % (indent, self.suffix(i, l), values))
        return lines

    def term_values(self, indent, vtype, jamp):
        """The bilinears of the jamps used in the quadratic forms (jamp(i) returns the C++ expression for jamp i)"""
        lines = []
        for (i, l) in self.terms:
            if i == l:
                expr = 'cxabs2( %s )' % jamp(i)
            else:
                expr = 'cxreal( %s ) * cxreal( %s ) + cximag( %s ) * cximag( %s )' % (jamp(i), jamp(l), jamp(i), jamp(l))
            lines.append('%sconst %s bsfA%s = %s;' % (indent, vtype, self.suffix(i, l), expr))
        return lines

    def link_sum(self):
        """The right hand side of the accumulation of one link"""
        parts = []
        for (i, l) in self.terms:
            s = self.suffix(i, l)
            if i == l: parts.append('cBsfS%s[ibsf] * bsfA%s' % (s, s))
            else: parts.append('2 * cBsfS%s[ibsf] * bsfA%s' % (s, s))
        return ' + '.join(parts)

    def ans2_arrays(self, indent):
        """The static constexpr arrays for ANS(2): the partner helicities and the colour matrix"""
        lines = []
        lines.append('%sstatic constexpr short cAns2Partner[nans2][ncomb] = {' % indent)
        for partner in self.ans2partner:
            lines.append('%s  { %s },' % (indent, ', '.join(['%d' % p for p in partner])))
        lines.append('%s};' % indent)
        lines.append('%sstatic constexpr fptype cAns2Cf[ncolor][ncolor] = {' % indent)
        for row in self.cf:
            lines.append('%s  { %s },' % (indent, ', '.join([format_fraction(Fraction(x)) for x in row])))
        lines.append('%s};' % indent)
        return lines

    def get_ans2_decl_cpp(self):
        """Before the C++ helicity loop of sigmaKin: the jamps of all good helicities of this page are kept for ANS(2)"""
        if self.nans2 == 0: return ''
        lines = []
        lines.append('      // AV vfold step 14g: the jamps of all good helicities of this page, for ANS(2)')
        lines.append('      cxtype_sv jampAll_sv[ncomb * nParity * ncolor];')
        return '\n'.join(lines) + '\n'

    def get_ans2_sum_cpp(self):
        """After the C++ helicity loop of sigmaKin: ANS(2) = sum_hel sum_ij CF(i,j) JAMP(+)_j conj(JAMP(-)_i)"""
        if self.nans2 == 0: return ''
        lines = []
        lines.append('      // AV vfold step 14g: ANS(2) of the Fortran SBORN, for each possible massless mother (see nans2)')
        lines.append('      {')
        lines += self.ans2_arrays('        ')
        lines.append('        int ighelOfHel[ncomb];')
        lines.append('        for( int ihel = 0; ihel < ncomb; ihel++ ) ighelOfHel[ihel] = -1;')
        lines.append('        for( int ighel = 0; ighel < cNGoodHel; ighel++ ) ighelOfHel[cGoodHel[ighel]] = ighel;')
        lines.append('        for( int iParity = 0; iParity < nParity; ++iParity )')
        lines.append('        {')
        lines.append('          fptype_sv* bsfx_sv = s_bsfsAccum + ( ( ievt00 + iParity * neppV ) / neppV ) * nbx;')
        lines.append('          for( int ians2 = 0; ians2 < nans2; ians2++ )')
        lines.append('          {')
        lines.append('            for( int ighel = 0; ighel < cNGoodHel; ighel++ )')
        lines.append('            {')
        lines.append('              const int ihel2 = cAns2Partner[ians2][cGoodHel[ighel]];')
        lines.append('              if( ihel2 < 0 ) continue; // not the BACK_HEL of this leg')
        lines.append('              const int ighel2 = ighelOfHel[ihel2];')
        lines.append('              if( ighel2 < 0 ) continue; // the other helicity vanishes')
        lines.append('              const cxtype_sv* jamp1_sv = jampAll_sv + ( ighel * nParity + iParity ) * ncolor;')
        lines.append('              const cxtype_sv* jamp2_sv = jampAll_sv + ( ighel2 * nParity + iParity ) * ncolor;')
        lines.append('              for( int icol = 0; icol < ncolor; icol++ )')
        lines.append('              {')
        lines.append('                for( int jcol = 0; jcol < ncolor; jcol++ )')
        lines.append('                {')
        lines.append('                  bsfx_sv[ibxAns2 + 2 * ians2] += cAns2Cf[icol][jcol] * ( cxreal( jamp2_sv[jcol] ) * cxreal( jamp1_sv[icol] ) + cximag( jamp2_sv[jcol] ) * cximag( jamp1_sv[icol] ) );')
        lines.append('                  bsfx_sv[ibxAns2 + 2 * ians2 + 1] += cAns2Cf[icol][jcol] * ( cximag( jamp2_sv[jcol] ) * cxreal( jamp1_sv[icol] ) - cxreal( jamp2_sv[jcol] ) * cximag( jamp1_sv[icol] ) );')
        lines.append('                }')
        lines.append('              }')
        lines.append('            }')
        lines.append('          }')
        lines.append('        }')
        lines.append('      }')
        return '\n'.join(lines) + '\n'

    #----------------------------------------------------------------------------------

    def get_definitions(self):
        """The definitions after cGoodHel in CPPProcess.cc"""
        ngroups = len(self.diagram_groups)
        lines = []
        lines.append('  // AV vfold step 9i+9j, %s only: the %s COLOUR-LINKED BORNS' % (self.process_string, number_word(self.nbsf)))
        lines.append('  // B_mn = <M|T_m.T_n|M> for the coloured legs m,n in {%s}, i.e. what the' % ','.join(['%d' % n for n in self.colored_legs]))
        lines.append('  // MadFKS Fortran routines B_SF_001..%03d compute from the cached Born amplitudes.' % self.nbsf)
        lines.append('  // Each B_SF_00N maps the %d amplitudes linearly onto JAMP1 (the plain Born colour' % self.namp)
        lines.append('  // flows, identical to jamp_sv here) and JAMP2 (the colour-linked flows), and')
        lines.append('  // returns Re(sum_ij CF(j,i) JAMP2(j) conj(JAMP1(i)))/IDEN summed over helicities.')
        lines.append('  // For this process JAMP2 is itself a fixed linear map of JAMP1 (both depend only')
        lines.append('  // on the %s diagram group%s %s), so each link collapses to a' % (number_word(ngroups), 's' if ngroups != 1 else '', join_words([format_ranges(g) for g in self.diagram_groups])))
        lines.append('  // quadratic form in jamp_sv with small %s coefficients (exact, derived with' % ('INTEGER' if self.all_integer else 'RATIONAL'))
        lines.append('  // rational arithmetic from the Fortran coefficient tables):')
        lines.append('  //   B_k = %s' % ' + '.join([self.term_formula(i, l) for (i, l) in self.terms]))
        lines.append('  // i.e. the same shape as the colour sum for |M|^2 (cf = %s), with a' % self.cf_string())
        lines.append('  // different matrix per link. Link k follows the Fortran order B_SF_00k, i.e. the')
        lines.append('  // upper triangle %s.' % ' '.join(['(%d,%d)' % tuple(link) for link in self.links]))
        lines.append('  // Evaluated per helicity right after the jamps, hence vectorized over the events of')
        lines.append('  // a SIMD page for free, and normalised exactly like the MEs (/helcolDenominators,')
        lines.append('  // but WITHOUT the multichannel numerator/denominator factor).')
        lines.append('  // [381 port: the new cudacpp computes the jamps in calculate_jamps and the colour sum')
        lines.append('  // in a separate color_sum function, so the colour-linked Borns and JAMP2 are now')
        lines.append('  // computed in sigmaKin from the jamps of each helicity, while AMP2 stays in calculate_jamps]')
        lines.append('  constexpr int nbsf = %d;' % self.nbsf)
        lines.append('  // AV vfold step 9i+9j: the same per-page buffer also carries,')
        lines.append('  // for each event, everything else that the Fortran SBORN returns and that')
        lines.append('  // its many callers need, so that mode C no longer has to run the scalar')
        lines.append('  // SBORN at all: the per-diagram |amp|^2 (Fortran AMP2, used by the')
        lines.append('  // multichannel single-diagram enhancement) and the per-colour-flow')
        lines.append('  // |jamp|^2 (Fortran JAMP2, used by the MC counterterms), both summed over')
        lines.append('  // helicities and NOT divided by the helicity/colour denominators (exactly')
        lines.append('  // like AMP2 and JAMP2 in the Fortran born.f).')
        lines.append('  constexpr int nbxDiag = CPPProcess::nbornamps; // %d (amplitudes, as in the Fortran AMP2)' % self.namp)
        lines.append('  // AV vfold step 14g: and the second component ANS(2) of the Fortran SBORN (real and imaginary parts), needed')
        lines.append('  // if the mother of the FKS splitting is massless, for the legs {%s} (normalised like the MEs).' % ','.join(['%d' % n for n in self.ans2legs]))
        lines.append('  constexpr int nans2 = %d;' % self.nans2)
        lines.append('  constexpr int ibxAns2 = nbsf + nbxDiag + ncolor;')
        lines.append('  constexpr int nbx = ibxAns2 + 2 * nans2; // per-event block: [bsf(%d), amp2(%d), jamp2(%d), ans2(2*%d)]' % (self.nbsf, self.namp, self.ncolor, self.nans2))
        lines.append('  static_assert( nbx == CPPProcess::nbornextras, "nbx must match CPPProcess::nbornextras" );')
        lines.append('  // Step 10a: also on GPUs. There the per-event block lives in a device buffer')
        lines.append('  // [nevt*nbx] owned by MatrixElementKernelDevice and passed to sigmaKin (argument allBsfx);')
        lines.append('  // it is registered after each launch, so that fbridgegetbsf_ and fbridgegetbornextras_')
        lines.append('  // can copy the values of the last call back to the host.')
        lines.append('  // [381 port: the new cudacpp runs one GPU stream per good helicity, so AMP2 is written by')
        lines.append('  // calculate_jamps into a per-helicity buffer ghelAllAmp2[nGoodHel][nbornamps][nevt] and the')
        lines.append('  // whole block is then summed over helicities in a fixed order in the compute_bornextras')
        lines.append('  // kernel, from ghelAllAmp2 and ghelAllJamps (no atomicAdd: results are reproducible)]')
        lines.append('#ifndef MGONGPUCPP_GPUIMPL')
        lines.append('  static fptype_sv* s_bsfs = nullptr;      // [npagV][nbx]: results of the last sigmaKin call (aligned raw buffer')
        lines.append('                                           // as in MemoryBuffers.h: std::vector<fptype_sv> ignores the alignment attribute)')
        lines.append('  static int s_bsfsSize = 0;               // npagV * nbx for the allocated s_bsfs')
        lines.append('  static fptype_sv* s_bsfsAccum = nullptr; // = s_bsfs inside sigmaKin, else nullptr (e.g. sigmaKin_getGoodHel)')
        lines.append('#else')
        lines.append('  static const fptype* s_devBsfs = nullptr; // device buffer [nevt*nbx] used by the last sigmaKin launch (host variable)')
        lines.append('#endif')
        lines.append('  static int s_bsfsNevt = 0; // nevt of the last sigmaKin call')
        return '\n'.join(lines) + '\n'

    def get_sigmakin_cpp(self):
        """The colour-linked Borns and JAMP2 in the C++ helicity loop of sigmaKin"""
        lines = []
        lines.append('        // *** COLOUR-LINKED BORNS BELOW *** (AV vfold, see nbsf above)')
        lines.append('        for( int iParity = 0; iParity < nParity; ++iParity )')
        lines.append('        {')
        if self.nbsf > 0:
            lines.append("          // [NB keep 'static' for these constexpr arrays, as for cf (issue #283)]")
            lines += self.coefficient_arrays('          ')
        lines.append('          const cxtype_sv* jampP_sv = jamp_sv + iParity * ncolor;')
        lines.append('          fptype_sv* bsfx_sv = s_bsfsAccum + ( ( ievt00 + iParity * neppV ) / neppV ) * nbx;')
        lines += self.term_values('          ', 'fptype_sv', lambda i: 'jampP_sv[%d]' % i)
        if self.nbsf > 0:
            lines.append('          for( int ibsf = 0; ibsf < nbsf; ibsf++ )')
            lines.append('            bsfx_sv[ibsf] += %s;' % self.link_sum())
        for i in range(self.ncolor):
            s = self.suffix(i, i)
            lines.append('          bsfx_sv[nbsf + nbxDiag + %d] += bsfA%s; // AV vfold: Fortran JAMP2(%d)' % (i, s, i + 1))
        if self.nans2 > 0:
            lines.append('          for( int icol = 0; icol < ncolor; icol++ ) jampAll_sv[( ighel * nParity + iParity ) * ncolor + icol] = jampP_sv[icol]; // for ANS(2)')
        lines.append('        }')
        return '\n'.join(lines) + '\n'

    def get_kernel_gpu(self):
        """The CUDA kernel summing the colour-linked Borns, AMP2 and JAMP2 over helicities (before sigmaKin)"""
        lines = []
        lines.append('#ifdef MGONGPUCPP_GPUIMPL')
        lines.append('  // AV vfold (step 10a): per-event block [bsf(%d), amp2(%d), jamp2(%d), ans2(2*%d)] summed over the good helicities' % (self.nbsf, self.namp, self.ncolor, self.nans2))
        lines.append('  // (in the order of cGoodHel, as in C++), from the jamps and per-diagram |amp|^2 of each helicity')
        lines.append('  __global__ void')
        lines.append('  compute_bornextras( fptype* allBsfx,              // output: per-event block [nevt*nbx]')
        lines.append('                      const fptype* ghelAllJamps,   // input: jamps super-buffer[2][ncolor][nGoodHel][nevt]')
        lines.append('                      const fptype* ghelAllAmp2,    // input: amp2 super-buffer[nGoodHel][nbornamps][nevt]')
        lines.append('                      const int helcolDenominator ) // input: helicity/colour denominator')
        lines.append('  {')
        lines.append('    using J_ACCESS = DeviceAccessJamp;')
        if self.nbsf > 0:
            lines.append("    // [NB keep 'static' for these constexpr arrays, as for cf (issue #283): needed in CUDA device code]")
            lines += self.coefficient_arrays('    ')
        if self.nans2 > 0:
            lines += self.ans2_arrays('    ')
        lines.append('    const int nevt = gridDim.x * blockDim.x;')
        lines.append('    const int ievt = blockDim.x * blockIdx.x + threadIdx.x;')
        lines.append('    fptype bsfx[nbx] = { 0 };')
        lines.append('    for( int ighel = 0; ighel < dcNGoodHel; ighel++ )')
        lines.append('    {')
        for i in range(self.ncolor):
            lines.append('      const cxtype jamp%d = J_ACCESS::kernelAccessIcolIhelNhelConst( ghelAllJamps, %d, ighel, dcNGoodHel );' % (i, i))
        lines.append('      for( int iamp = 0; iamp < nbxDiag; iamp++ )')
        lines.append('        bsfx[nbsf + iamp] += ghelAllAmp2[( (size_t)ighel * nbxDiag + iamp ) * nevt + ievt]; // Fortran AMP2')
        lines += self.term_values('      ', 'fptype', lambda i: 'jamp%d' % i)
        if self.nbsf > 0:
            lines.append('      for( int ibsf = 0; ibsf < nbsf; ibsf++ )')
            lines.append('        bsfx[ibsf] += %s;' % self.link_sum())
        for i in range(self.ncolor):
            s = self.suffix(i, i)
            lines.append('      bsfx[nbsf + nbxDiag + %d] += bsfA%s; // Fortran JAMP2(%d)' % (i, s, i + 1))
        lines.append('    }')
        if self.nans2 > 0:
            lines.append('    // Step 14g: ANS(2) of the Fortran SBORN, for each possible massless mother (see nans2)')
            lines.append('    int ighelOfHel[ncomb];')
            lines.append('    for( int ihel = 0; ihel < ncomb; ihel++ ) ighelOfHel[ihel] = -1;')
            lines.append('    for( int ighel = 0; ighel < dcNGoodHel; ighel++ ) ighelOfHel[dcGoodHel[ighel]] = ighel;')
            lines.append('    for( int ians2 = 0; ians2 < nans2; ians2++ )')
            lines.append('    {')
            lines.append('      for( int ighel = 0; ighel < dcNGoodHel; ighel++ )')
            lines.append('      {')
            lines.append('        const int ihel2 = cAns2Partner[ians2][dcGoodHel[ighel]];')
            lines.append('        if( ihel2 < 0 ) continue; // not the BACK_HEL of this leg')
            lines.append('        const int ighel2 = ighelOfHel[ihel2];')
            lines.append('        if( ighel2 < 0 ) continue; // the other helicity vanishes')
            lines.append('        for( int icol = 0; icol < ncolor; icol++ )')
            lines.append('        {')
            lines.append('          const cxtype jamp1 = J_ACCESS::kernelAccessIcolIhelNhelConst( ghelAllJamps, icol, ighel, dcNGoodHel );')
            lines.append('          for( int jcol = 0; jcol < ncolor; jcol++ )')
            lines.append('          {')
            lines.append('            const cxtype jamp2 = J_ACCESS::kernelAccessIcolIhelNhelConst( ghelAllJamps, jcol, ighel2, dcNGoodHel );')
            lines.append('            bsfx[ibxAns2 + 2 * ians2] += cAns2Cf[icol][jcol] * ( cxreal( jamp2 ) * cxreal( jamp1 ) + cximag( jamp2 ) * cximag( jamp1 ) );')
            lines.append('            bsfx[ibxAns2 + 2 * ians2 + 1] += cAns2Cf[icol][jcol] * ( cximag( jamp2 ) * cxreal( jamp1 ) - cxreal( jamp2 ) * cximag( jamp1 ) );')
            lines.append('          }')
            lines.append('        }')
            lines.append('      }')
            lines.append('    }')
        lines.append('    // Final colour-linked Borns (same denominators as the MEs, no multichannel factor; amp2 and jamp2 are raw sums)')
        lines.append('    for( int ibsf = 0; ibsf < nbsf; ibsf++ ) bsfx[ibsf] /= helcolDenominator;')
        lines.append('    for( int ibx = ibxAns2; ibx < nbx; ibx++ ) bsfx[ibx] /= helcolDenominator; // ANS(2)')
        lines.append('    for( int ibx = 0; ibx < nbx; ibx++ ) allBsfx[ievt * nbx + ibx] = bsfx[ibx];')
        lines.append('  }')
        lines.append('#endif')
        lines.append('')
        lines.append('  //--------------------------------------------------------------------------')
        lines.append('')
        return '\n'.join(lines) + '\n'

    def get_fortran_interface(self):
        """The extern C functions for Fortran at the end of CPPProcess.cc"""
        lines = []
        lines.append('')
        lines.append('//==========================================================================')
        lines.append('')
        lines.append('// AV vfold step 9i+9j, %s only: Fortran access to the %s' % (self.process_string, number_word(self.nbsf)))
        lines.append('// colour-linked Borns of the events of the LAST sequence call on this library.')
        lines.append('// CALL FBRIDGEGETBSF(NEVT, BSFS, STATUS) with DOUBLE PRECISION BSFS(%d,NEVT):' % self.nbsf)
        lines.append('// STATUS=0 on success; 1 if NEVT does not match the last sigmaKin call.')
        lines.append('// Step 10a: exported by the CUDA library too (copy from the device).')
        lines.append('#ifdef MGONGPUCPP_GPUIMPL')
        lines.append('namespace mg5amcBsf = mg5amcGpu;')
        lines.append('#else')
        lines.append('namespace mg5amcBsf = mg5amcCpu;')
        lines.append('#endif')
        lines.append('extern "C"')
        lines.append('{')
        lines.append('  void fbridgegetbsf_( const int* pnevt, double* bsfs, int* pstatus )')
        lines.append('  {')
        lines += self.status_check()
        lines.append('    {')
        lines.append('      *pstatus = 1;')
        lines.append('      return;')
        lines.append('    }')
        lines.append('    mg5amcBsf::getColourLinkedBorns( *pnevt, bsfs );')
        lines.append('    *pstatus = 0;')
        lines.append('  }')
        lines.append('')
        lines.append('  // AV vfold (step 9i+9j): the whole per-event block of the LAST sequence call:')
        lines.append('  // CALL FBRIDGEGETBORNEXTRAS(NEVT, EXTRAS, STATUS), DOUBLE PRECISION EXTRAS(%d,NEVT)' % self.nbx)
        lines.append('  // = the %d colour-linked Borns, then AMP2(1:%d), then JAMP2(1:%d), then ANS(2) (see nbx).' % (self.nbsf, self.namp, self.ncolor))
        lines.append('  // Step 10a: exported by the CUDA library too (copy from the device).')
        lines.append('  void fbridgegetbornextras_( const int* pnevt, double* extras, int* pstatus )')
        lines.append('  {')
        lines += self.status_check()
        lines.append('    {')
        lines.append('      *pstatus = 1;')
        lines.append('      return;')
        lines.append('    }')
        lines.append('    mg5amcBsf::getBornExtras( *pnevt, extras );')
        lines.append('    *pstatus = 0;')
        lines.append('  }')
        lines.append('')
        lines.append('  // AV vfold (step 14g): the size of the block returned by fbridgegetbornextras_, so that')
        lines.append('  // the Fortran side can check that it was generated with the same version of this code')
        lines.append('  void fbridgegetnbornextras_( int* pnbx )')
        lines.append('  {')
        lines.append('    *pnbx = mg5amcBsf::nbx;')
        lines.append('  }')
        lines.append('}')
        lines.append('')
        return '\n'.join(lines)

    @staticmethod
    def status_check():
        return ['#ifdef MGONGPUCPP_GPUIMPL',
                '    if( *pnevt != mg5amcBsf::s_bsfsNevt || mg5amcBsf::s_devBsfs == nullptr )',
                '#else',
                '    if( *pnevt != mg5amcBsf::s_bsfsNevt )',
                '#endif']
