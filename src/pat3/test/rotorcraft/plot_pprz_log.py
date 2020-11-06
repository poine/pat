#! /usr/bin/env python

import logging, sys, os, time, math, numpy as np
import matplotlib.pyplot as plt
import pdb

import pat3.pprz_bag as p3_ppbag

# <message name="ROTORCRAFT_FP" id="147">
#       <field name="east"     type="int32" alt_unit="m" alt_unit_coef="0.0039063"/>
#       <field name="north"    type="int32" alt_unit="m" alt_unit_coef="0.0039063"/>
#       <field name="up"       type="int32" alt_unit="m" alt_unit_coef="0.0039063"/>
#       <field name="veast"    type="int32" alt_unit="m/s" alt_unit_coef="0.0000019"/>
#       <field name="vnorth"   type="int32" alt_unit="m/s" alt_unit_coef="0.0000019"/>
#       <field name="vup"      type="int32" alt_unit="m/s" alt_unit_coef="0.0000019"/>
#       <field name="phi"      type="int32" alt_unit="deg" alt_unit_coef="0.0139882"/>
#       <field name="theta"    type="int32" alt_unit="deg" alt_unit_coef="0.0139882"/>
#       <field name="psi"      type="int32" alt_unit="deg" alt_unit_coef="0.0139882"/>
#       <field name="carrot_east"   type="int32" alt_unit="m" alt_unit_coef="0.0039063"/>
#       <field name="carrot_north"  type="int32" alt_unit="m" alt_unit_coef="0.0039063"/>
#       <field name="carrot_up"     type="int32" alt_unit="m" alt_unit_coef="0.0039063"/>
#       <field name="carrot_psi"    type="int32" alt_unit="deg" alt_unit_coef="0.0139882"/>
#       <field name="thrust"        type="int32"/>
#       <field name="flight_time"   type="uint16" unit="s"/>
# </message>

class MyDataSet(p3_ppbag.DataSet):
    def __init__(self):
        def rec_gps_int(msg): return [float(msg.ecef_x)/100, float(msg.ecef_y)/100, float(msg.ecef_z)/100]
        def rec_rotorcraft_fp(msg):
            return [float(msg.east)*0.0039063, float(msg.north)*0.0039063, float(msg.up)*0.0039063]
        p3_ppbag.DataSet.__init__(self, msg_ids=['ROTORCRAFT_FP'], msg_cbks=[rec_rotorcraft_fp])
        

        

def main(filename):
    d = MyDataSet()
    d.readPprzBag(filename)
    plt.plot(d.tss[0], d.vals[0])
    plt.show()
    
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    filename = '/home/poine/work/pat/data/logs/20_08_04__13_41_08.data'
    main(filename)
