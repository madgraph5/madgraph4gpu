# Copyright (C) 2023-2024 CERN.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: Z. Wettersten (Sep 2024) for the MG5aMC CUDACPP plugin.

import os
import subprocess
import re
import sys
import importlib.util
SPEC_EXPORTCPP = importlib.util.find_spec('madgraph.iolibs.export_cpp')
PLUGIN_export_cpp = importlib.util.module_from_spec(SPEC_EXPORTCPP)
SPEC_EXPORTCPP.loader.exec_module(PLUGIN_export_cpp)
sys.modules['PLUGIN.CUDACPP_OUTPUT.PLUGIN_export_cpp'] = PLUGIN_export_cpp # allow 'import PLUGIN.CUDACPP_OUTPUT.PLUGIN_export_cpp' in model_handling.py
del SPEC_EXPORTCPP
###print('id(export_cpp)=%s'%id(export_cpp))
###print('id(PLUGIN_export_cpp)=%s'%id(PLUGIN_export_cpp))

# AV - use template files from PLUGINDIR instead of MG5DIR
###from madgraph import MG5DIR
PLUGINDIR = os.path.dirname( __file__ )

# AV - model_handling includes the custom FileWriter, ALOHAWriter, UFOModelConverter, OneProcessExporter and HelasCallWriter, plus additional patches
import PLUGIN.CUDACPP_OUTPUT.model_handling as model_handling
import PLUGIN.CUDACPP_OUTPUT.output as output

# AV - create a plugin-specific logger
import logging
logger = logging.getLogger('madgraph.PLUGIN.CUDACPP_OUTPUT.output')
from madgraph import MG5DIR
#------------------------------------------------------------------------------------

from os.path import join as pjoin
import madgraph
import madgraph.iolibs.files as files
import madgraph.iolibs.export_v4 as export_v4
import madgraph.various.misc as misc
import madgraph.interface.reweight_interface as rwgt_interface
import madgraph.various.banner as banner
import models.check_param_card as check_param_card
import madgraph.interface.extended_cmd as extended_cmd
import madgraph.interface.common_run_interface as common_run_interface

from . import launch_plugin

class TREX_OneProcessExporter(model_handling.PLUGIN_OneProcessExporter):
    """A custom OneProcessExporter for the TREX reweighting"""
    
    rwgt_template = 'gpu/rwgt_runner.inc'

    # ZW - rwgt functions
    def get_rwgt_legs(self, process):
        """Return string with particle ids and status in the REX std::pair format"""
        return ",".join(["{\"%i\",\"%i\"}" % (leg.get('state'), leg.get('id')) \
            for leg in process.get('legs')]).replace('0', '-1')
    
    def get_rwgt_legs_vec(self, processes):
        """Return string with vectors of particle ids and statuses"""
        prtSets = []
        for k in range(len(processes)):
            prtSets.append("{" + self.get_rwgt_legs(processes[k]) + "}")
        return ",".join(prtSets)
        
    def get_init_prts_vec(self, process):
        """Return string with initial state particle ids for use in REX event sorting"""
        prts = ",".join(["\"%i\"" % leg.get('id') for leg in process.get('legs') if leg.get('state') == 0])
        return "{" + prts + "}"
    
    def get_init_prts_vecs(self, processes):
        """Return string with vectors of initial state particle ids"""
        prtSets = []
        for k in range(len(processes)):
            prtSets.append(self.get_init_prts_vec(processes[k]))
        return ",".join(prtSets)
    
    def get_fin_prts_vec(self, process):
        """Return string with final state particle ids for use in REX event sorting"""
        prts = ",".join(["\"%i\"" % leg.get('id') for leg in process.get('legs') if leg.get('state') == 1])
        return "{" + prts + "}"
    
    def get_fin_prts_vecs(self, processes):
        """Return string with vectors of final state particle ids"""
        prtSets = []
        for k in range(len(processes)):
            prtSets.append(self.get_fin_prts_vec(processes[k]))
        return ",".join(prtSets)
        
    def get_rwgt_procMap(self, process):
        """Return string with particle states and order in the REX procMap format"""
        currState = False
        retString = "thisProc{{\"-1\",{"
        for leg in process.get('legs'):
            if currState == leg.get('state'):
                retString += "\"%i\"," % leg.get('id')
            else:
                currState = leg.get('state')
                retString += "}},{\"1\",{\"%i,\"" % leg.get('id')
        retString = retString[:-1] + "}}}"
        return retString
    
    def get_proc_dir(self):
        """Return process directory name for the current process"""
        return "P%d_%s" % (self.process_number, self.process_name)
        
    def get_rwgt_runner(self):
        """Return string to initialise the rwgtRunners in teawREX"""
        return "%s::runner" % (self.get_proc_dir())
    
    def get_rwgt_includes(self):
        """Return string with the include directives for the REX reweighting"""
        return "#include \"P%d_%s/rwgt_runner.cc\"" % (self.process_number, self.process_name)
    
    def write_rwgt_header(self):
        """Writes a simple rwgt_runner.h file to forward declare the runner object"""
        # Adjust the placeholders for use with `.format()`
        rwgt_h = """#ifndef {namespace}_RWGT_RUNNER_H
    #define {namespace}_RWGT_RUNNER_H
    #include \"rwgt_instance.h\"
    namespace {namespace} {{
        extern rwgt::instance runner;
    }}
    #endif""".format(namespace=self.get_proc_dir())
        
        # Using `with` statement for better file handling
        with open(os.path.join(self.path, 'rwgt_runner.h'), 'w') as ff:
            ff.write(rwgt_h)
    
    def edit_rwgt_header(self):
        """Adds process-specific details to the rwgt_runner.h template"""
        replace_dict = super().get_process_class_definitions(write=False)
        replace_dict['process_namespace'] = self.get_proc_dir()
        replace_dict['info_lines'] = model_handling.PLUGIN_export_cpp.get_mg5_info_lines()
        template = open(pjoin(self.template_path,'REX', 'rwgt_runner.h'),'r').read()
        ff = open(pjoin(self.path, 'rwgt_runner.h'),'w')
        ff.write(template % replace_dict)
        ff.close()
    
    def edit_rwgt_runner(self):
        """Create the rwgt_runner.cc file for the REX reweighting"""
        ###misc.sprint('Entering PLUGIN_OneProcessExporterRwgt.edit_rwgt_runner')
        # Create the rwgt_runner.cc file
#        replace_dict = {}
        replace_dict = super().get_process_class_definitions(write=False)
#        rwgt_runner = self.get_proc_dir() + self.rwgt_template
        replace_dict['process_namespace'] = self.get_proc_dir()
        replace_dict['info_lines'] = model_handling.PLUGIN_export_cpp.get_mg5_info_lines()
        replace_dict['init_prt_ids'] = self.get_init_prts_vecs(self.matrix_elements[0].get('processes'))
        replace_dict['fin_prt_ids'] = self.get_fin_prts_vecs(self.matrix_elements[0].get('processes'))
        replace_dict['process_events'] = self.get_rwgt_legs_vec(self.matrix_elements[0].get('processes'))
        replace_dict['no_events'] = len(self.matrix_elements[0].get('processes'))
        template = open(pjoin(self.template_path,'REX', 'rwgt_runner.inc'),'r').read()
        ff = open(pjoin(self.path, 'rwgt_runner.cc'),'w')
        ff.write(template % replace_dict)
        ff.close()
    
    # ZW - override the PLUGIN method to generate the rwgt_runner.cc file as well
    # note: also generating standard check_sa.cc and gcheck_sa.cu files, which
    # are not used in the REX reweighting
    def generate_process_files(self):
        """Generate mgOnGpuConfig.h, CPPProcess.cc, CPPProcess.h, check_sa.cc, gXXX.cu links"""
        # misc.sprint('Entering RWGT_OneProcessExporter.generate_process_files')
        super().generate_process_files()
        # misc.sprint('Generating rwgt_runner files')
        self.edit_rwgt_header()
        self.edit_rwgt_runner()
        # misc.sprint('Finished generating rwgt files')
        
class TREX_ProcessExporter(output.PLUGIN_ProcessExporter):
    
    oneprocessclass = TREX_OneProcessExporter
    
    rwgt_names = []
    proc_lines = []
    
    s = PLUGINDIR + '/madgraph/iolibs/template_files/'
    from_template = dict(output.PLUGIN_ProcessExporter.from_template)
    from_template['src'] = from_template['src'] + [s+'REX/REX.cc', s+'REX/teawREX.cc',
                                                   s+'REX/REX.h', s+'REX/teawREX.h',
                                                    s+'REX/rwgt_instance.h', s+'REX/rwgt_instance.cc']
    from_template['SubProcesses'] = from_template['SubProcesses'] + [s+'gpu/cudacpp_driver.mk',
                                                                     s+'REX/rwgt_instance.h', s+'REX/REX.h', s+'REX/teawREX.h']

    to_link_in_P = output.PLUGIN_ProcessExporter.to_link_in_P + ['rwgt_instance.h', 'REX.h', 'teawREX.h']

    template_src_make = pjoin(PLUGINDIR, 'madgraph' ,'iolibs', 'template_files','gpu','cudacpp_rex_src.mk')
    template_tst_make = pjoin(PLUGINDIR, 'madgraph', 'iolibs', 'template_files','gpu','cudacpp_test.mk')
    template_Sub_make = pjoin(PLUGINDIR, 'madgraph', 'iolibs', 'template_files','gpu','cudacpp_runner.mk')
    
    def generate_subprocess_directory(self, matrix_element, cpp_helas_call_writer,
                                      proc_number=None):
        """Generate the Pxxxxx directory for a subprocess in C++ standalone,
        including the necessary .h and .cc files"""

        
        process_exporter_cpp = self.oneprocessclass(matrix_element,cpp_helas_call_writer)
        
        self.rwgt_names.append("P%d_%s" % (process_exporter_cpp.process_number, 
                                             process_exporter_cpp.process_name))
        
        process_lines = "\n".join([process_exporter_cpp.get_process_info_lines(me) for me in \
                                   process_exporter_cpp.matrix_elements])
        self.proc_lines.append(process_lines)
        
        # Create the directory PN_xx_xxxxx in the specified path
        dirpath = pjoin(self.dir_path, 'SubProcesses', "P%d_%s" % (process_exporter_cpp.process_number, 
                                             process_exporter_cpp.process_name))
        try:
            os.mkdir(dirpath)
        except os.error as error:
            logger.warning(error.strerror + " " + dirpath)
    
        with misc.chdir(dirpath):
            logger.info('Creating files in directory %s' % dirpath)
            process_exporter_cpp.path = dirpath
            # Create the process .h and .cc files
            process_exporter_cpp.generate_process_files()
            for file in self.to_link_in_P:
                files.ln('../%s' % file) 
        return
    
    def export_driver(self):
        # misc.sprint("In export_driver")
        # misc.sprint("Current working directory is: %s" % self.dir_path)
        replace_dict = {}
        replace_dict['info_lines'] = model_handling.PLUGIN_export_cpp.get_mg5_info_lines()
        replace_dict['multiprocess_lines'] = "\n".join(self.proc_lines)
        replace_dict['include_lines'] = ''
        replace_dict['run_set'] = ''
        replace_dict['fbridge_vec'] = ''
        for name in self.rwgt_names:
            replace_dict['include_lines'] += '#include "%s/rwgt_runner.h"\n' % name
            replace_dict['run_set'] += '%s::getEventSet(),' % name
            replace_dict['fbridge_vec'] += '%s::bridgeConstr(),' % name
        replace_dict['run_set'] = replace_dict['run_set'][:-1]
        replace_dict['fbridge_vec'] = replace_dict['fbridge_vec'][:-1]
        template_path = os.path.join( PLUGINDIR, 'madgraph', 'iolibs', 'template_files' )
        template = open(pjoin(template_path,'REX', 'rwgt_driver.inc'),'r').read()
        ff = open(pjoin(self.dir_path, 'SubProcesses', 'rwgt_driver.cc'),'w')
        ff.write(template % replace_dict)
        ff.close()
    
    def link_makefile(self):
        """Link the makefile for the REX reweighting"""
        files.ln(pjoin(self.dir_path, 'SubProcesses', 'cudacpp_driver.mk'), starting_dir=pjoin(self.dir_path, 'SubProcesses'), name='makefile')
    
    def finalize(self, matrix_element, cmdhistory, MG5options, outputflag):
        self.export_driver()
        self.link_makefile()
        return super().finalize(matrix_element, cmdhistory, MG5options, outputflag)
    
class TREX_ReweightInterface(rwgt_interface.ReweightInterface):
    """A custom ReweightInterface for the TREX reweighting"""
    
    sa_class = 'standalone_trex'
    
    def __init__(self, *args, **kwargs):
        """Initialise the TREX reweighting interface
        Currently no (substantial) changes compared to upstream are necessary,
        but adding an __init__ method allows for future modifications"""
        super().__init__(*args, **kwargs)
        self.debug_output = 'tREX_debug'
        self.param_card = None
        self.reweight_card = []
        self.reweight_names = []
        
    def setup_f2py_interface(self):
        """"Override native setup_f2py_interface to avoid parsing things not necessary for TREX reweighting"""
        
        self.create_standalone_directory()
        self.compile()
        
    def launch_actual_reweighting(self, *args, **kwargs):
        """override standard launch command to instead call the TREX reweighting"""
        
        import csv
        
        if self.rwgt_dir:
            path_me =self.rwgt_dir
        else:
            path_me = self.me_dir 
        
        if self.second_model or self.second_process or self.dedicated_path:
            rw_dir = pjoin(path_me, 'rw_me_%s' % self.nb_library)
        else:
            rw_dir = pjoin(path_me, 'rw_me')

        run_path = pjoin(rw_dir, 'SubProcesses')
        input_file = os.path.relpath(self.lhe_input.path, run_path)
        output_file = input_file + 'rw'
        output_path = self.lhe_input.path + 'rw'
        param_card = pjoin(rw_dir, 'Cards', 'param_card.dat')

        #ZW: Exceptions, making sure all the necessary files for teawREX are accessible
        if( misc.is_executable(pjoin(run_path,'rwgt_driver_gpu.exe')) ):
            driver = pjoin(run_path, 'rwgt_driver_gpu.exe')
        elif(misc.is_executable(pjoin(run_path,'rwgt_driver_cpp.exe')) ):
            driver = pjoin(run_path,'rwgt_driver_cpp.exe')
        else:
            raise Exception('No teawREX driver found for parallel reweighting')
        if not os.path.exists(param_card):
            try:
                files.cp(os.path.join(path_me, 'Cards', 'param_card_default.dat'), param_card)
            except:
                raise Exception("No param_card.dat file found in %s" % pjoin(path_me, 'Cards'))
        param_path = os.path.relpath(param_card, run_path)    

        rwgt_card = os.path.join(path_me, 'Cards', 'reweight_card.dat')
        
        self.write_reweight_card(rwgt_card)

        if not os.path.exists(rwgt_card):
            try:
                files.cp(os.path.join(path_me, 'Cards', 'reweight_card_default.dat'), rwgt_card)
            except:
                raise Exception("No reweight_card.dat file found in %s" % pjoin(path_me, 'Cards'))
        rwgt_path = os.path.relpath(rwgt_card, run_path)
        target = ''
        if not self.mother:
            name, ext = self.lhe_input.name.rsplit('.',1)
            target = '%s_out.%s' % (name, ext)            
        elif self.output_type != "default" :
            target = pjoin(self.mother.me_dir, 'Events', self.mother.run_name, 'events.lhe')
        else:
            target = self.lhe_input.path

        #ZW: rwgt_driver is written and compiled properly, now just to figure out how to run it through MG
        subprocess.call([driver, '-lhe=%s' % input_file, '-slha=%s' % param_card, '-rwgt=%s' % rwgt_card, '-out=%s' % output_file], cwd=run_path)

        files.mv(output_path, target)
        csv_file = pjoin(run_path, 'rwgt_results.csv')
        with open(csv_file, newline='') as results:
            iters = csv.reader(results)
            for row in iters:
                self.all_cross_section[(row[0],'')] = (float(row[1]), float(row[2]))

        return
    
    def compile(self):
        """override compile to use the TREX makefiles"""
        
        if self.multicore=='wait':
            return
        
        if not self.rwgt_dir:
            path_me = self.me_dir
        else:
            path_me = self.rwgt_dir
        
        rwgt_dir_possibility =   ['rw_me','rw_me_%s' % self.nb_library,'rw_mevirt','rw_mevirt_%s' % self.nb_library]
        for onedir in rwgt_dir_possibility:
            if not os.path.isdir(pjoin(path_me,onedir)):
                continue
            pdir = pjoin(path_me, onedir, 'SubProcesses')
            if self.mother:
                nb_core = self.mother.options['nb_core'] if self.mother.options['run_mode'] !=0 else 1
            else:
                nb_core = 1
            misc.compile(cwd=pdir, nb_core=nb_core,mode='cpp')
        return
    
    def load_module(self):
        """override load_module since we do not use it"""
        return
    
    # def import_command_file(self, filepath):
    #     """override import_command_file to simply launch TREX"""
    #     self.exec_cmd('launch', precmd=True)
    #     return
    
    def do_launch(self, line):
        """override do_launch to instead overwrite the reweight_card
        to fit the expected input for TREX without having to extend TREX itself"""
        args = self.split_arg(line)
        opts = self.check_launch(args)
        mgcmd = self.mg5cmd
        if opts['rwgt_name']:
            self.options['rwgt_name'] = opts['rwgt_name']
        if opts['rwgt_info']:
            self.options['rwgt_info'] = opts['rwgt_info']
        model_line = self.banner.get('proc_card', 'full_model_line')

        # TV: Load model: needed for the combine_ij function: maybe not needed everyt time??? 
        model = self.banner.get('proc_card', 'model')
        self.load_model( model, True, False)

        if not self.has_standalone_dir:
            out = self.setup_f2py_interface()
            if out:
                return

        if not self.param_card:
            s_orig = self.banner['slha']
            self.param_card = check_param_card.ParamCard(s_orig.splitlines())

        # get the mode of reweighting #LO/NLO/NLO_tree/...
        type_rwgt = self.get_weight_names()
        
        if self.rwgt_dir:
            path_me =self.rwgt_dir
        else:
            path_me = self.me_dir 
        
            
        # get iterator over param_card and the name associated to the current reweighting.
        param_card_iterator, tag_name = self.handle_param_card(model_line, args, type_rwgt)

        self.reweight_names.append(tag_name)
        
        # perform the scanning
        if param_card_iterator:
            if self.options['rwgt_name']:
                reweight_name = self.options['rwgt_name'].rsplit('_',1)[0] # to avoid side effect during the scan
            else:
                reweight_name = None
            for i,card in enumerate(param_card_iterator):
                if reweight_name:
                    self.options['rwgt_name'] = '%s_%s' % (reweight_name, i+1)
                self.new_param_card = card
                #card.write(pjoin(rw_dir, 'Cards', 'param_card.dat'))
                self.exec_cmd("launch --keep_card", printcmd=False, precmd=True)
    
    def check_multicore(self):
        """override check_multicore to overloading the CPU (we never want to run TREX in multicore mode)"""
        return False
    
    def handle_param_card(self, model_line, args, type_rwgt):
        """override handle_param_card to get rid of all the unnecessary checks and file writing
        now simply loads the param_card and uses get_diff to tranlate into internal format"""    

        if self.rwgt_dir:
            path_me =self.rwgt_dir
        else:
            path_me = self.me_dir 
            
        if self.second_model or self.second_process or self.dedicated_path:
            rw_dir = pjoin(path_me, 'rw_me_%s' % self.nb_library)
        else:
            rw_dir = pjoin(path_me, 'rw_me')
        if not '--keep_card' in args:
            if self.has_nlo and self.rwgt_mode != "LO":
                rwdir_virt = rw_dir.replace('rw_me', 'rw_mevirt')
            with open(pjoin(rw_dir, 'Cards', 'param_card.dat'), 'w') as fsock:
                fsock.write(self.banner['slha']) 
            out, cmd = common_run_interface.CommonRunCmd.ask_edit_card_static(cards=['param_card.dat'],
                                ask=self.ask, pwd=rw_dir, first_cmd=self.stored_line,
                                write_file=False, return_instance=True
                                )
            self.stored_line = None
            card = cmd.param_card
            new_card = card.write()
        elif self.new_param_card:
            new_card = self.new_param_card.write()
        else:
            new_card = open(pjoin(rw_dir, 'Cards', 'param_card.dat')).read()
        
        # check for potential scan in the new card 
        pattern_scan = re.compile(r'''^(decay)?[\s\d]*scan''', re.I+re.M) 
        param_card_iterator = []
        if pattern_scan.search(new_card):
            import madgraph.interface.extended_cmd as extended_cmd
            try:
                import internal.extended_cmd as extended_internal
                Shell_internal = extended_internal.CmdShell
            except:
                Shell_internal = None
            if not isinstance(self.mother, (extended_cmd.CmdShell, Shell_internal)): 
                raise Exception("scan are not allowed on the Web")
            # at least one scan parameter found. create an iterator to go trough the cards
            main_card = check_param_card.ParamCardIterator(new_card)
            if self.options['rwgt_name']:
                self.options['rwgt_name'] = '%s_0' % self.options['rwgt_name']

            param_card_iterator = main_card
            first_card = param_card_iterator.next(autostart=True)
            new_card = first_card.write()
            self.new_param_card = first_card
            #first_card.write(pjoin(rw_dir, 'Cards', 'param_card.dat'))  

        # check if "Auto" is present for a width parameter)
        if 'block' not in new_card.lower():
            raise Exception(str(new_card))
        tmp_card = new_card.lower().split('block',1)[1]
        if "auto" in tmp_card: 
            if param_card_iterator:
                first_card.write(pjoin(rw_dir, 'Cards', 'param_card.dat'))
            else:
                ff = open(pjoin(rw_dir, 'Cards', 'param_card.dat'),'w')
                ff.write(new_card)
                ff.close()
                
            self.mother.check_param_card(pjoin(rw_dir, 'Cards', 'param_card.dat'))
            new_card = open(pjoin(rw_dir, 'Cards', 'param_card.dat')).read()


        # Find new tag in the banner and add information if needed
        if 'initrwgt' in self.banner and self.output_type == 'default': 
            if 'name=\'mg_reweighting\'' in self.banner['initrwgt']:
                blockpat = re.compile(r'''<weightgroup name=\'mg_reweighting\'\s*weight_name_strategy=\'includeIdInWeightName\'>(?P<text>.*?)</weightgroup>''', re.I+re.M+re.S)
                before, content, after = blockpat.split(self.banner['initrwgt'])
                header_rwgt_other = before + after
                pattern = re.compile('<weight id=\'(?:rwgt_(?P<id>\\d+)|(?P<id2>[_\\w\\-\\.]+))(?P<rwgttype>\\s*|_\\w+)\'>(?P<info>.*?)</weight>', re.S+re.I+re.M)
                mg_rwgt_info = pattern.findall(content)
                maxid = 0
                for k,(i, fulltag, nlotype, diff) in enumerate(mg_rwgt_info):
                    if i:
                        if int(i) > maxid:
                            maxid = int(i)
                        mg_rwgt_info[k] = (i, nlotype, diff) # remove the pointless fulltag tag
                    else:
                        mg_rwgt_info[k] = (fulltag, nlotype, diff) # remove the pointless id tag
                        
                maxid += 1
                rewgtid = maxid
                if self.options['rwgt_name']:
                    #ensure that the entry is not already define if so overwrites it
                    for (i, nlotype, diff) in mg_rwgt_info[:]:
                        for flag in type_rwgt:
                            if 'rwgt_%s' % i == '%s%s' %(self.options['rwgt_name'],flag) or \
                                i == '%s%s' % (self.options['rwgt_name'], flag):
                                    logger.warning("tag %s%s already defines, will replace it", self.options['rwgt_name'],flag)
                                    mg_rwgt_info.remove((i, nlotype, diff))
                                                
            else:
                header_rwgt_other = self.banner['initrwgt'] 
                mg_rwgt_info = []
                rewgtid = 1
        else:
            self.banner['initrwgt']  = ''
            header_rwgt_other = ''
            mg_rwgt_info = []
            rewgtid = 1

        # add the reweighting in the banner information:
        #starts by computing the difference in the cards.
        #s_orig = self.banner['slha']
        #self.orig_param_card_text = s_orig
        s_new = new_card
        self.new_param_card = check_param_card.ParamCard(s_new.splitlines())
        
        #define tag for the run
        if self.options['rwgt_name']:
            tag = self.options['rwgt_name']
        else:
            tag = str(rewgtid)

        if 'rwgt_info' in self.options and self.options['rwgt_info']:
            card_diff = self.options['rwgt_info']
            for name in type_rwgt:
                mg_rwgt_info.append((tag, name, self.options['rwgt_info']))
        elif not self.second_model and not self.dedicated_path:
            old_param = self.param_card
            new_param =  self.new_param_card
            card_diff = old_param.create_diff(new_param)
            if card_diff == '' and not self.second_process:
                    logger.warning(' REWEIGHTING: original card and new card are identical.')
            try:
                if old_param['sminputs'].get(3)- new_param['sminputs'].get(3) > 1e-3 * new_param['sminputs'].get(3):
                    logger.warning("We found different value of alpha_s. Note that the value of alpha_s used is the one associate with the event and not the one from the cards.")
            except Exception as error:
                logger.debug("error in check of alphas: %s" % str(error))
                pass #this is a security                
            if not self.second_process:
                for name in type_rwgt:
                    mg_rwgt_info.append((tag, name, card_diff))
            else:
                str_proc = "\n change process  ".join([""]+self.second_process)
                for name in type_rwgt:
                    mg_rwgt_info.append((tag, name, str_proc + '\n'+ card_diff))
        else:
            if self.second_model:
                str_info = "change model %s" % self.second_model
            else:
                str_info =''
            if self.second_process:
                str_info += "\n change process  ".join([""]+self.second_process)
            if self.dedicated_path:
                for k,v in self.dedicated_path.items():
                    str_info += "\n change %s %s" % (k,v)
            card_diff = str_info
            str_info += '\n' + s_new
            for name in type_rwgt:
                mg_rwgt_info.append((tag, name, str_info))

        # re-create the banner.
        self.banner['initrwgt'] = header_rwgt_other
        if self.output_type == 'default':
            self.banner['initrwgt'] += '\n<weightgroup name=\'mg_reweighting\' weight_name_strategy=\'includeIdInWeightName\'>\n'
        else:
            self.banner['initrwgt'] += '\n<weightgroup name=\'main\'>\n'
        for tag, rwgttype, diff in mg_rwgt_info:
            if self.inc_sudakov:
                try:
                    sud_order = int(rwgttype[-1]) -1
                    sud_order = '10' +rwgttype[-2:]
                    self.banner['initrwgt'] += '<weight id=\'%s\'>%sscale_%s_sud</weight>\n' % \
                            (rwgttype, diff, sud_order)
                except IndexError:
                    logger.critical('This is a reweighted event file! Do not reweight with ewsudakov twice')
                    sys.exit(1)
            else:
                if tag.isdigit():
                    self.banner['initrwgt'] += '<weight id=\'rwgt_%s%s\'>%s</weight>\n' % \
                                    (tag, rwgttype, diff)
                else:
                    self.banner['initrwgt'] += '<weight id=\'%s%s\'>%s</weight>\n' % \
                                    (tag, rwgttype, diff)
        self.banner['initrwgt'] += '\n</weightgroup>\n'
        self.banner['initrwgt'] = self.banner['initrwgt'].replace('\n\n', '\n')

        #logger.info('starts to compute weight for events with the following modification to the param_card:')
        #logger.info(card_diff.replace('\n','\nKEEP:'))
        try:
            self.run_card = banner.Banner(self.banner).charge_card('run_card')
        except Exception:
            logger.debug('no run card found -- reweight interface')
            self.run_card = None

        if self.options['rwgt_name']:
            tag_name = self.options['rwgt_name']
        else:
            tag_name = 'rwgt_%s' % rewgtid

        self.reweight_card.append(card_diff)
        
        return param_card_iterator, tag_name
    
    def write_reweight_card(self,rwgt_path):
        """function for collecting all the reweight iterations from the parsed reweight card
        and write it out with the explicit 'set BLOCK PARAM VALUE' format"""
        if( len(self.reweight_names) != len(self.reweight_card) ):
            raise Exception('Mismatch in number of reweight names and reweight cards')
        
        output_card = ''
        
        for i, card in enumerate(self.reweight_card):
            output_card += 'launch --rwgt_name=%s\n' % self.reweight_names[i]
            output_card += card + '\n'
        
        output_card = output_card.replace('param_card', '').replace('  ', ' ')
        
        with open(rwgt_path, 'w') as f:
            f.write(output_card)
        
        return
    
    def do_quit(self, line):
        if self.exitted:
            return
        
        self.launch_actual_reweighting()
        
        self.exitted = True
        
        if 'init' in self.banner:
            cross = 0 
            error = 0
            for line in self.banner['init'].split('\n'):
                split = line.split()
                if len(split) == 4:
                    cross += float(split[0])
                    error += float(split[1])**2
            error = error**0.5
        if not self.multicore == 'create':
            # No print of results for the multicore mode for the one printed on screen
            if 'orig' not in self.all_cross_section:
                logger.info('Original cross-section: %s +- %s pb' % (cross, error))
            else: 
                logger.info('Original cross-section: %s +- %s pb (cross-section from sum of weights: %s)' % (cross, error, self.all_cross_section['orig'][0]))
            logger.info('Computed cross-section:')
            keys = list(self.all_cross_section.keys())
            keys.sort(key=lambda x: str(x))
            for key in keys:
                if key == 'orig':
                    continue
                logger.info('%s : %s +- %s pb' % (key[0] if not key[1] else '%s%s' % key,
                    self.all_cross_section[key][0],self.all_cross_section[key][1] ))  
        self.terminate_fortran_executables()

        if self.rwgt_dir and self.multicore == False:
            self.save_to_pickle()
        
        with misc.stdchannel_redirected(sys.stdout, os.devnull):
            for run_id in self.calculator:
                del self.calculator[run_id]
            del self.calculator