#! /usr/bin/env python
#-*- coding: utf-8 -*-

import numpy as np
import scipy.integrate
from io import StringIO
import xml.etree.ElementTree as ET
import logging, ast

"""
 This is a 6dof model for a fixed wing vehicle
"""
LOG = logging.getLogger(__name__)
import pdb


class Param:
    """
    Represents the parameters of the flight dynamic model
    """
    def __init__(self, filename=None, avl_data=None):
        self.g = 9.81
        if filename is not None : self.read_from_xml(filename)
        elif avl_data is not None : self.compute_from_avl(avl_data) # TODO...
        self.compute_auxiliary()

    def compute_auxiliary(self):
        """
        compute auxiliary coefficients
        """
        self.invJ = np.linalg.inv(self.J)
        self.I, self.invI = self.J, self.invJ # FIXME
        try:
            self.inv_mamat = np.linalg.inv(self.mamat)
        except AttributeError:
            pass
        self.eng_to_body=[]
        for i in range(0, self.eng_nb):
            self.eng_to_body.append(np.eye(3))
        self.input_nb = self.eng_nb+self.sfc_nb
        for c in ["CL", "CY", "Cl", "Cm", "Cn"]:
            val = [getattr(self, '{:s}_{:s}'.format(c, ax)) for ax in ['p', 'q', 'r']]
            setattr(self, c+'_omega', val)

    def read_from_xml(self, filename):
        """
          Quick and dirty xml parsing
        """
        LOG.info("\n  Loading parameters from {:s}".format(filename))
        tree = ET.parse(filename)
        root = tree.getroot()

        self.name = root.attrib["name"]

        masses = root.find("masses")
        self.m = float(masses.find("mass").text)
        try:
            self.J = np.genfromtxt(StringIO(masses.find("inertia").text), delimiter=",")
        except TypeError:   # python 2 compatibility
            self.J = np.genfromtxt(StringIO(unicode(masses.find("inertia").text)), delimiter=",")
            
        aerodynamics = root.find("aerodynamics")
        for param in ["Sref", "Cref", "Bref"]:
            setattr(self, param, float(aerodynamics.find(param).text))

        stability = aerodynamics.find("stability")
        misc = stability.find("misc")
        for att in ['Vref', 'CL0', 'alpha0', 'CD0', 'CD_k1', 'CD_k2', 'Cm0']:
            setattr(self, att, float(misc.attrib[att]))
        for s in ['alpha', 'beta', 'p', 'q', 'r']:
            x = stability.find(s)
            for c in ["CL", "CY", "Cl", "Cm", "Cn"]:
                setattr(self, c+'_'+s, float(x.attrib[c]))
            
        surfaces = stability.findall("surface");
        self.sfc_name = [s.attrib["name"] for s in surfaces]
        self.sfc_nb = len(self.sfc_name)
        for c in ["CL", "CY", "CD", "Cl", "Cm", "Cn"]:
            setattr(self, c+"_sfc", np.array([float(s.attrib[c]) for s in surfaces]))

        propulsion = root.find("propulsion")
        engines = propulsion.findall("engine")
        self.eng_name = [e.attrib["name"] for e in engines]
        self.eng_nb = len(self.eng_name)
        try:
            setattr(self, "eng_pos", [np.genfromtxt(StringIO(e.attrib["pos"]),delimiter=",") for e in engines])
            setattr(self, "eng_align", [np.genfromtxt(StringIO(e.attrib["align"]),delimiter=",") for e in engines])
        except TypeError:   # python 2 compatibility
            setattr(self, "eng_pos", [np.genfromtxt(StringIO(unicode(e.attrib["pos"])),delimiter=",") for e in engines])
            setattr(self, "eng_align", [np.genfromtxt(StringIO(unicode(e.attrib["align"])),delimiter=",") for e in engines])
        for param in ["fmax", "rhoi", "nrho", "Vi", "nV", "tau"]:
            setattr(self, param+"s", [float(e.attrib[param]) for e in engines])
        
        try:
            control = root.find("control")
            trim_q = control.find('trim_q')
            code = trim_q.attrib["code"]
            self.trim_Usfc_of_elevator = code
            #https://greentreesnakes.readthedocs.org/en/latest/index.html
            #https://www.ianlewis.org/en/dynamically-adding-method-classes-or-class-instanc
            #pdb.set_trace()
        except AttributeError:
            self.trim_Usfc_of_elevator = None


    def compute_from_avl(self, avl_data):
        self.name = avl_data['name']

        self.m = avl_data['inertial']['mass']
        self.J = avl_data['inertial']['inertia']

        v_CL0 = avl_data['coef_CL0']
        v_CD0 = avl_data['coef_CD0']
        v_ref  = avl_data['coef_ref']
        sfc = avl_data['surfaces']

        for param in ['Sref', 'Cref', 'Bref']:
            setattr(self, param, float(v_ref[param]))
        self.CL0 = float(v_CL0['CLtot']); self.alpha0 = 0.
        self.CD0 = float(v_CD0['CDtot'])
        self.CD_k1 = 0.
        self.CD_k2 = (float(v_ref['CDtot'])-float(v_CD0['CDtot']))/float(v_ref['CLtot'])**2
        #K2 = 1/(pi*Bref**2/Sref*e)
        self.Cm0 = float(v_CL0['Cmtot'])
        for ax in ['alpha', 'beta', 'p', 'q', 'r']:
            for coef in ['CL', 'CY', 'Cl', 'Cm', 'Cn']:
                setattr(self, '{:s}_{:s}'.format(coef, ax), float(v_ref['{:s}{:s}'.format(coef, ax[0])]))

        # reordering... we want aileron, elevator, rudder, flaps, then the rest
        pat_sfc = ['aileron','elevator', 'rudder', 'flap']
        my_sfc = list(sfc); self.sfc_nb = len(my_sfc)
        idxs = []
        for s in pat_sfc:
            try: idxs.append(my_sfc.index(s))
            except: pass
        for _i in range(self.sfc_nb):
            if _i not in idxs:
                idxs.append(_i)
                pat_sfc.append(my_sfc[_i])
        self.sfc_name = pat_sfc[0:self.sfc_nb]
        for coef, coef_avl in [('CL', 'CL'), ('CY','CY'), ('CD','CDff'), ('Cl','Cl'), ('Cm','Cm'), ('Cn', 'Cn')]:
            val = [np.rad2deg(float(v_ref['{:s}d{:d}'.format(coef_avl, idxs[i]+1)])) for i in range(0,self.sfc_nb)]
            setattr(self, '{:s}_sfc'.format(coef), val)

        prop = avl_data['propulsion']
        self.eng_nb = len(prop['name'])
        self.eng_name = prop['name']
        self.eng_pos = np.array(prop['pos']); self.eng_align = np.array(prop['align']); self.fmaxs = prop['fmax']
        self.rhois = prop['rhoi']; self.nrhos = prop['nrho']
        self.Vis = prop['Vi']; self.nVs = prop['nV']
        self.taus = prop['tau']
        self.Vref = 15. # FIXME

    def save_xml(self, filename):
        """
        Quick and dirty xml printing
        """
        LOG.info(" saving fdm configuration to {}".format(filename))
        with open(filename, 'w') as f:
            f.write('''<?xml version="1.0"?>
<fdm_config name="{:s}">\n'''.format(self.name))
            f.write('''
  <masses>
    <mass unit="kg"> {:f} </mass>
    <inertia unit="kg.m2">
'''.format(self.m))
            np.savetxt(f, self.I, '%.4f', delimiter=', ')
            f.write('''    </inertia>
  </masses>\n\n''')
            f.write('  <propulsion>\n')
            for i in range(0, self.eng_nb):
                pos_str = StringIO(); np.savetxt(pos_str, self.eng_pos[i][np.newaxis], '%.2f', delimiter=', ')
                align_str = StringIO(); np.savetxt(align_str, self.eng_align[i][np.newaxis], '%.2f', delimiter=', ')
                fmt = '    <engine name="{:s}" pos="{:s}" align="{:s}" fmax="{:f}" rhoi="{:f}" nrho="{:f}" Vi="{:f}" nV="{:f}" tau="{:f}"/>\n'
                f.write(fmt.format(self.eng_name[i], pos_str.getvalue().strip(), align_str.getvalue().strip(), self.fmaxs[i],
                                   self.rhois[i], self.nrhos[i], self.Vis[i], self.nVs[i], self.taus[i]))
            f.write('  </propulsion>\n\n')
            f.write('  <aerodynamics>\n')
            f.write('    <Sref unit="m2" comment="Reference surface">{:f}</Sref>\n'.format(self.Sref))
            f.write('    <Cref unit="m"  comment="Reference chord">  {:f}</Cref>\n'.format(self.Cref))
            f.write('    <Bref unit="m"  comment="Reference length"> {:f}</Bref>\n'.format(self.Bref))
            f.write('    <stability>\n')
            fmt = '       <misc Vref="{:f}" CL0="{:f}" alpha0="{:f}" CD0="{:f}" CD_k1="{:f}" CD_k2="{:f}" Cm0="{:f}"/>\n'
            f.write(fmt.format(self.Vref, self.CL0, self.alpha0, self.CD0, self.CD_k1, self.CD_k2, self.Cm0))
            fmt = '       <alpha                   CL="{: f}" CY="{: f}"                Cl="{: f}" Cm="{: f}" Cn="{: f}"/>\n'
            f.write(fmt.format(self.CL_alpha, self.CY_alpha, self.Cl_alpha, self.Cm_alpha, self.Cn_alpha))
            fmt = '       <beta                    CL="{: f}" CY="{: f}"                Cl="{: f}" Cm="{: f}" Cn="{: f}"/>\n'
            f.write(fmt.format(self.CL_beta, self.CY_beta, self.Cl_beta, self.Cm_beta, self.Cn_beta))
            fmt = '       <p                       CL="{: f}" CY="{: f}"                Cl="{: f}" Cm="{: f}" Cn="{: f}"/>\n'
            f.write(fmt.format(self.CL_p, self.CY_p, self.Cl_p, self.Cm_p, self.Cn_p))
            fmt = '       <q                       CL="{: f}" CY="{: f}"                Cl="{: f}" Cm="{: f}" Cn="{: f}"/>\n'
            f.write(fmt.format(self.CL_q, self.CY_q, self.Cl_q, self.Cm_q, self.Cn_q))
            fmt = '       <r                       CL="{: f}" CY="{: f}"                Cl="{: f}" Cm="{: f}" Cn="{: f}"/>\n'
            f.write(fmt.format(self.CL_r, self.CY_r, self.Cl_r, self.Cm_r, self.Cn_r))
            fmt = '       <surface name="{:s}"{:s} CL="{: f}" CY="{: f}" CD="{: f}" Cl="{: f}" Cm="{: f}" Cn="{: f}"/>\n'
            for i in range(0,self.sfc_nb):
                f.write(fmt.format(self.sfc_name[i], "".ljust(8-len(self.sfc_name[i])),
                                   self.CL_sfc[i], self.CY_sfc[i], self.CD_sfc[i],
                                   self.Cl_sfc[i], self.Cm_sfc[i], self.Cn_sfc[i]))
            f.write('    </stability>\n')
            f.write('  </aerodynamics>\n')
            f.write('</fdm_config>\n')


        
    def __str__(self):
        """
          Quick and dirty human readable printing of parameters
        """
        res  = "Name: {:s}".format(self.name)
        res += """
  Mass: {:f} kg
  Inertia:\n {:s} kg.m2
        """.format(self.m, str(self.J))
        res += """
  Aero:
    Sref: {:f} m2    
    Cref: {:f} m    
    Bref: {:f} m    
    """.format(self.Sref, self.Cref, self.Bref)
        for m in ['Vref', 'CL0', 'alpha0', 'CD0', 'CD_k1', 'CD_k2', 'Cm0']: res += (m+": {:f} ".format(getattr(self, m)))
        res+="\n           alpha      beta         p         q         r\n"
        for c in ["CL", "CY", "Cl", "Cm", "Cn"]:
            res+= "    {:s}".format(c)
            for s in ['alpha', 'beta', 'p', 'q', 'r']:
                res+=" {: f}".format(getattr(self, c+"_"+s))
            res+="\n"
        res += "      "
        for name in self.sfc_name: res+= "{:>10s}".format(name)
        res+="\n"
        for c in ["CL", "CY", "CD", "Cl", "Cm", "Cn"]:
            res+= "    {:s}".format(c)
            for i in range(0, len(self.sfc_name)):
                res+=" {: f}".format(getattr(self, c+"_sfc")[i])
            res+="\n"

        res += "\n  Propulsion:\n"
        res += "          "
        for name in self.eng_name: res+= "{:>16s}".format(name)
        res+="\n"

        for param in ["pos", "align"]:
            res+="    {:8s}".format(param)
            for i in range(0, self.eng_nb):
                #pdb.set_trace()
                #res += "{:>16s}".format(getattr(self, "eng_"+param)[i])
                res += "{:>16s}".format(str(getattr(self, "eng_"+param)[i]))
            res+="\n"

        for param in ["fmax", "rhoi", "nrho", "Vi", "nV", "tau"]:
            res += "    {:8s}".format(param)
            for i in range(0, self.eng_nb):
                res += "{:16f}".format(getattr(self, param+"s")[i])
            res+="\n"
        return res

    # some default actuator allocation
    def u_slice_eng(self): return slice(0, self.eng_nb)
    def u_slice_sfc(self): return slice(self.eng_nb, self.eng_nb+self.sfc_nb)
    def u_size(self): return self.input_nb
    def u_elevator(self): return self.eng_nb+1  # let's assume elevator is the second surface
