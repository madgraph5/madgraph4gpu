################################################################################
#
# Copyright (c) 2009 The MadGraph5_aMC@NLO Development team and Contributors
#
# This file is a part of the MadGraph5_aMC@NLO project, an application which 
# automatically generates Feynman diagrams and matrix elements for arbitrary
# high-energy processes in the Standard Model and beyond.
#
# It is subject to the MadGraph5_aMC@NLO license which should accompany this 
# distribution.
#
# For more information, visit madgraph.phys.ucl.ac.be and amcatnlo.web.cern.ch
#
################################################################################
from __future__ import division
from __future__ import absolute_import
import subprocess
import unittest
import os
import re
import shutil
import sys
import logging
import time
import tempfile
import math
import madgraph


logger = logging.getLogger('test_cmd')

import tests.unit_tests.iolibs.test_file_writers as test_file_writers

import madgraph.interface.master_interface as MGCmd
import madgraph.interface.madevent_interface as MECmd
import madgraph.interface.launch_ext_program as launch_ext
import madgraph.iolibs.files as files

import madgraph.various.misc as misc
import madgraph.various.lhe_parser as lhe_parser
import madgraph.various.banner as banner_mod
import madgraph.various.lhe_parser as lhe_parser
import madgraph.various.banner as banner

_file_path = os.path.split(os.path.dirname(os.path.realpath(__file__)))[0]
_pickle_path =os.path.join(_file_path, 'input_files')

from madgraph import MG4DIR, MG5DIR, MadGraph5Error, InvalidCmd

from tests.acceptance_tests.test_cmd_madevent import *

pjoin = os.path.join


#===============================================================================
# TestCmd
#===============================================================================
class TestCPPfromfile(TestMEfromfile): # inherit from upstream test_cmd_madevent
    """test that we can launch everything from a single file"""

 

    def test_ggtt_mixed(self):
        """checking time of flight is working fine"""

        if logging.getLogger('madgraph').level <= 20:
            stdout=None
            stderr=None
        else:
            devnull =open(os.devnull,'w')
            stdout=devnull
            stderr=devnull
            
        try:
            shutil.rmtree('/tmp/MGPROCESS/')
        except Exception as error:
            pass
        
        cmd = """import model sm
                 set automatic_html_opening False --no_save
                 set notification_center False --no_save
                 generate p p > t t~
                 output madevent_simd %s -f -nojpeg
                 launch  
                 set nevents 100
                 """ %self.run_dir

        open(pjoin(self.path, 'mg5_cmd'),'w').write(cmd)
        
        subprocess.call([sys.executable, pjoin(MG5DIR, 'bin','mg5_aMC'), 
                         pjoin(self.path, 'mg5_cmd')],
                         #cwd=self.path,
                        stdout=stdout, stderr=stderr)

        self.check_parton_output(cross=505.5, error=2.749)
        event = '%s/Events/run_01/unweighted_events.lhe' % self.run_dir
        if not os.path.exists(event):
            misc.gunzip(event)
        
        has_zero = False
        has_non_zero = False
        lhefile = lhe_parser.EventFile(event)
        lhefile.apply_fct_on_event(fcts=lhe_parser.Event.check) 

        nb_event = 0
        for event in lhe_parser.EventFile(event):
            event.check()
            nb_event+=1

        self.assertEqual(nb_event, 100)
        
        self.assertFalse(self.debuging)
    

