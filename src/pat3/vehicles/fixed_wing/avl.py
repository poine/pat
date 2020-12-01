#!/usr/bin/env python
#-*- coding: utf-8 -*-
#
# Copyright 2013-2020 Antoine Drouin (poinix@gmail.com)
#
# This file is part of PAT.
#
#    PAT is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    PAT is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with PAT.  If not, see <http://www.gnu.org/licenses/>.
#

"""
Support for Avl (Athena Vortex Lattice)
See: http://web.mit.edu/drela/Public/web/avl/

This program is used to compute aerodynamic coefficients (and incidentally inertias)  
for our fixed wing flight dynamic model by using Avl (digital windtunnel...).

1: How does it work
This version works by running avl binary interactive text interface (i.e., the regular avl program) in a child subprocess.
It means, if you are able to run 'avl' and interract with it from your terminal, this
program should work.

I have cython based version that allows to call internal avl (fortran) routines from python,
but this requires avl sources for building the interface, which makes it less convenient,
albeit more flexible, robust and performant. Feel free to contact me if interested.

2: How does it work
We run 3 experiments on Avl, two for finding zero lift and zero drag conditions, and one at a 'reference speed'

3: Does it work?
Kindoff... I am partly convinced that what we do here is not completely incorrect.
Extracting the correct coefficients from Avl requires an exact understanding of Avl assumptions and internals.
I have a bunch of (broken) unit test from long ago. Again, feel free to contact me if interested.

"""
import os, sys, logging, pexpect, re, math, numpy as np, ast
import pdb

LOG = logging.getLogger(__name__)

#import pat2.solid as ps, pat2.misc_utils as mu
#import pat2.vehicle.fixed_wing as pv
import pat3.vehicles.fixed_wing.legacy_6dof as p3_v

s_u, s_w, s_q, s_the, s_v, s_p, s_r, s_phi, s_x, s_y, s_z, s_psi, s_size = np.arange(0,13)
avl2s1 = [s_x, s_y, s_z, s_u, s_v, s_w, s_phi, s_the, s_psi, s_p, s_q, s_r]
#tex_s = ps.tex_s1[avl2s1]

def get_run_cases(filename):
    f, cases, case = open(filename), [], None
    for line in f:
        m=re.match(' Run case\s+(\S+):\s+(\S+)\s+', line)
        if m:
            num, name = m.groups()
            case = {'name':name}
            cases.append(case)
        else:
            for _re, name in [[' velocity  =\s+(\S+)\s+m/s', 'velocity'],
                             [' alpha     =\s+(\S+)\s+deg', 'alpha']]:
                m=re.match(_re, line)
                if m: case[name] = float(m.groups()[0])
    return cases

def get_surfaces(filename):
    f, control_surfaces, found_control = open(filename), [], False
    for line in f:
        if line.strip().startswith('CONTROL'): found_control = True
        elif found_control:
            ctl_name = line.strip().split(' ')[0]
            if not ctl_name in control_surfaces: control_surfaces.append(ctl_name)
            found_control = False

    return control_surfaces


def change_to_infile_dir_and_update_outfile(infile, outfile=None):
    infile_dir = os.path.dirname(infile)            # avl needs to run from the config files directory
    org_dir = os.getcwd()                           # save that to be able to return
    cfg_file = os.path.basename(infile)             # strip the directory part of input file
    if outfile != None:
        abspath_outfile = os.path.abspath(outfile)  # get absolute path to output file
        if os.path.isfile(abspath_outfile):         # avl won't erase its outfile if it exists,
            os.unlink(abspath_outfile)              # so we do it
    os.chdir(infile_dir)                            #
    if outfile!=None: res = cfg_file, org_dir, abspath_outfile
    else: res = cfg_file, org_dir
    return res

def chat_with_avl(cfg_file, dialog, debug=False):
    child = pexpect.spawn ('avl '+cfg_file)
    if debug: child.logfile = sys.stdout
    for expect, reply in dialog:
        child.expect(expect)
        if reply is not None: child.sendline (reply)
    if debug: print(str(child))

def run_st(infile, outfile, case, debug=False):
    LOG.info(" Running ST  on {} runcase {} to {}".format(infile, case, outfile))
    cfg_file, org_dir, updated_outfile = change_to_infile_dir_and_update_outfile(infile, outfile)
    dialog = [['AVL   c>',                                           'OPER'],
              ['OPER \(case \S+\)   c>',                             '{:d}'.format(case)],
              ['OPER \(case \S+\)   c>',                             'X'],
              ['OPER \(case \S+\)   c>',                             'ST'],
              ['Enter filename, or <return> for screen output   s>', updated_outfile],
              ['OPER \(case \S+\)   c>',                             None]]
    chat_with_avl(cfg_file, dialog, debug)
    os.chdir(org_dir)

def run_eig(infile, outfile, case, debug=False):
    LOG.info(" Running EIG on %s runcase %d to %s", infile, case, outfile)
    cfg_file, org_dir, updated_outfile = change_to_infile_dir_and_update_outfile(infile, outfile)
    dialog = [['AVL   c>',                               'MODE'],
              ['MODE   c>',                              '{:d}'.format(case)],
              ['MODE   c>',                              'S'],
              ['Enter output filename \(or <Return>\):', updated_outfile],
              ['MODE   c>',                              None]]
    chat_with_avl(cfg_file, dialog, debug)
    os.chdir(org_dir)

# FIXME: this needs some love
def run_masses(infile, debug=False):
    """
    Runs Avl for parsing a mass file.
    returns computed mass and inertia
    """
    LOG.info(" run_masses on {}".format(infile))
    cfg_file, org_dir = change_to_infile_dir_and_update_outfile(infile)
    child = pexpect.spawn ('avl '+cfg_file)
    if debug: child.logfile = sys.stdout
    child.expect(r'Mass \s+=\s+(\S+)\s+kg')
    (mass,) = child.match.groups()
    def read_intertia():
        child.expect(r'Ixx -Ixy -Ixz   \|\s+(\S+)\s+(\S+)\s+(\S+)')
        (Ixx, mIxy, mIxz) = child.match.groups()
        child.expect(r'Iyy -Iyz = \|\s+(\S+)\s+(\S+)')
        (Iyy, mIyz) = child.match.groups()
        child.expect(r'Izz   \|\s+(\S+)')
        (Izz,) = child.match.groups()
        return Ixx, mIxy, mIxz, Iyy, mIyz, Izz
    
    Ixx, mIxy, mIxz, Iyy, mIyz, Izz = read_intertia()

    ###
    child.expect(r'mxx  mxy  mxz   \|\s+(\S+)\s+(\S+)\s+(\S+)')
    (mxx, mxy, mxz) = child.match.groups()
    child.expect(r'     myy  myz = \|\s+(\S+)\s+(\S+)')
    (myy, myz) = child.match.groups()
    child.expect(r'          mzz   \|\s+(\S+)')
    (mzz,) = child.match.groups()
    ###
    Iaxx, mIaxy, mIaxz, Iayy, mIayz, Iazz = read_intertia()

    child.expect ('AVL   c>')
    child.sendline ('Quit')
    child.expect(pexpect.EOF)
    I = np.array([[ float(Ixx),  float(mIxy), float(mIxz)],
                  [ float(mIxy), float(Iyy),  float(mIyz)],
                  [ float(mIxz), float(mIyz), float(Izz)]])
    Ia = np.array([[ float(Iaxx),  float(mIaxy), float(mIaxz)],
                  [ float(mIaxy), float(Iayy),  float(mIayz)],
                  [ float(mIaxz), float(mIayz), float(Iazz)]])
    mamat = np.array([[ float(mxx),  float(mxy), float(mxz)],
                      [ float(mxy), float(myy),  float(myz)],
                      [ float(mxz), float(myz), float(mzz)]])
    res = {'mass':float(mass), 'inertia':I, 'ap_inertia':Ia, 'mamat':mamat}
    os.chdir(org_dir)
    return res


def parse_st(filename, surfaces, debug=False):
    """
    Extracts coefficients from an avl st file
    """
    parser = [(r'.*Sref =\s+(\S+)\s+Cref =\s+(\S+)\s+Bref =\s+(\S+)',('Sref', 'Cref', 'Bref')),
              (r'.*Alpha =\s+(\S+)\s+', ('Alpha',)),
              (r'.*CYtot =\s+(\S+)\s+ Cmtot =\s+(\S+)\s+',('CYtot', 'Cmtot')),
              (r'.*CLtot =\s+(\S+)\s+',('CLtot',)),
              (r'.*CDtot =\s+(\S+)\s+',('CDtot',)),
              (r'.*CDvis =\s+(\S+)\s+ CDind =\s+(\S+)\s+',('CDvis', 'CDind')),
              (r'.*CLff  = \s+(\S+)\s+ CDff  =\s+(\S+)\s+',('CLff', 'CDff')),
              # surfaces will be inserted here...              
              (r'.*CLa =\s+(\S+)\s+ CLb =\s+(\S+)\s+', ('CLa', 'CLb')),
              (r'.*CYa =\s+(\S+)\s+ CYb =\s+(\S+)\s+', ('CYa', 'CYb')),
              (r'.*Cla =\s+(\S+)\s+ Clb =\s+(\S+)\s+', ('Cla', 'Clb')),
              (r'.*Cma =\s+(\S+)\s+ Cmb =\s+(\S+)\s+', ('Cma', 'Cmb')),
              (r'.*Cna =\s+(\S+)\s+ Cnb =\s+(\S+)\s+', ('Cna', 'Cnb')),
              (r'.*CLp =\s+(\S+)\s+ CLq =\s+(\S+)\s+ CLr =\s+(\S+)\s+', ('CLp', 'CLq', 'CLr')),
              (r'.*CYp =\s+(\S+)\s+ CYq =\s+(\S+)\s+ CYr =\s+(\S+)\s+', ('CYp', 'CYq', 'CYr')),
              (r'.*Clp =\s+(\S+)\s+ Clq =\s+(\S+)\s+ Clr =\s+(\S+)\s+', ('Clp', 'Clq', 'Clr')),
              (r'.*Cmp =\s+(\S+)\s+ Cmq =\s+(\S+)\s+ Cmr =\s+(\S+)\s+', ('Cmp', 'Cmq', 'Cmr')),
              (r'.*Cnp =\s+(\S+)\s+ Cnq =\s+(\S+)\s+ Cnr =\s+(\S+)\s+', ('Cnp', 'Cnq', 'Cnr')),
              # coef for surfaces will be appended here
              ]
    for sfc in reversed(surfaces):
        parser.insert(7, (r'.*{:s}\s+ =\s+(\S+)\s+'.format(sfc), (sfc,)))

    for coef in ['CL', 'CY', 'Cl', 'Cm', 'Cn', 'CDff', 'e']:
        regexp = ".+"+"".join(['{:s}d{:d} =\s+(\S+)\s+'.format(coef, i+1) for i,s in enumerate(surfaces)])
        var = tuple(['{:s}d{:d}'.format(coef, i+1) for i,s in enumerate(surfaces)])
        parser.append((regexp, var))
    
    f = open(filename)
    i = 0; values = {}
    for line in f:
        LOG.debug(" %s", line.replace('\n', ''))
        if i < len(parser):
            m=re.match(parser[i][0], line)
            if m:
                LOG.debug(" match found %s %s", parser[i][1], m.groups())
                for j, v in enumerate(m.groups()):
                    values[parser[i][1][j]] = v
                i += 1
    if i==len(parser): LOG.info(" successfully parsed {} ({} values)".format(filename, i))
    else: LOG.error(" failed to parse {} (failed after {} values on {})".format(filename, i, len(parser)))
    return values


def parse_eig(filename, va, alpha):
    ''' loads an AVL sysmat (aka jacobians) file and converts it to pat s1 format'''
    f = open(filename, 'r')
    f.readline();
    sv_txt, input_txt = f.readline().split('|')
    #print input_txt
    d = np.loadtxt(f)
    Aa, Ba = d[:, 0:s_size], d[:, s_size:]
    Ba = mu.deg_of_rad(Ba) # avl has control surfaces deflection in degres...
    A, B = Aa[:,avl2s1][avl2s1], Ba[avl2s1]
    #        x, y, z, u                   v  w                   phi theta psi p q r
    Xtrim = [0, 0, 0, va*math.cos(alpha), 0, va*math.sin(alpha), 0, alpha, 0, 0, 0, 0] 
    return A, B, Xtrim, Aa, Ba



'''
We assume that if foo.avl is the avl config file, we also have foo.run, foo.mass, foo.prop in the same directory
'''
def ac_name(avl_infile): return os.path.splitext(os.path.basename(avl_infile))[0]
def runfile(avl_infile): return os.path.splitext(avl_infile)[0]+'.run'
def massfile(avl_infile): return os.path.splitext(avl_infile)[0]+'.mass'
def propfile(avl_infile): return os.path.splitext(avl_infile)[0]+'.prop'
'''
AVL outputs: 
'''
def run_outfile(outdir, ac_name, run_name): return os.path.join(outdir, '{}_{}'.format(ac_name, run_name)) 
def st_outfile(outdir, ac_name, run_name):  return run_outfile(outdir, ac_name, run_name)+'.st'
def eig_outfile(outdir, ac_name, run_name): return run_outfile(outdir, ac_name, run_name)+'.eig'


def make_pat_config(avl_infile, outdir, _run_st=False, _run_eig=False, debug=False):
    LOG.info(" making PAT config from {}".format(avl_infile))
    _ac_name, _runfile = ac_name(avl_infile), runfile(avl_infile)
    runs = get_run_cases(_runfile)
    LOG.info("  found runcases {} in {}".format([c['name'] for c in runs], _runfile))
    for i, r in enumerate(runs):
        st_of = st_outfile(outdir, _ac_name, r['name'])
        if _run_st or not os.path.isfile(st_of): run_st(avl_infile, st_of,  i+1, debug)
        eig_of = eig_outfile(outdir, _ac_name, r['name'])
        if (_run_eig  or not os.path.isfile(eig_of)) and i>=2: run_eig(avl_infile, eig_of,  i+1, debug)

    surfaces = get_surfaces(avl_infile)
    val_CL0, val_CD0, val_ref = [parse_st(st_outfile(outdir, _ac_name, runs[i]['name']), surfaces, debug) for i in range(0, 3)]
    inertial = run_masses(avl_infile, debug=debug)
    f = open(propfile(avl_infile))
    propulsion = ast.literal_eval(f.read())
    avl_data = {'name': _ac_name,
                'inertial': inertial,
                'coef_CL0': val_CL0, 'coef_CD0': val_CD0, 'coef_ref': val_ref,
                'Vref': runs[2]['velocity'],
                'surfaces': surfaces,
                'propulsion': propulsion
                }
    return p3_v.Param(avl_data=avl_data)



