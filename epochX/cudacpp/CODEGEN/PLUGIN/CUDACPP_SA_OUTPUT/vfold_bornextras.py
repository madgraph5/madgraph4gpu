# Copyright (C) 2020-2025 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Oct 2026) for the MG5aMC CUDACPP plugin.
#
# AV road B (vfold): the --vfold output option, which enables the generation of the vfold extensions

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
