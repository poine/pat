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
        if filename != None : self.read_from_xml(filename)
        #elif avl_data <> None : self.compute_from_avl(avl_data) # TODO...
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
        LOG.info("loading parameters from {:s}".format(filename))
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
                res += "{:>16s}".format(getattr(self, "eng_"+param)[i])
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
