#! /usr/bin/env python
import logging, os, math, numpy as np, matplotlib.pyplot as plt
import pdb

#
# /mnt/mint17/home/poine/work/python-aerospace-toolbox/config/avl/cularis.avl
# 
# /mnt/mint17/home/poine/work/ngfw/trunk/simulations/py_avl/config/avl/cularis.avl
# ( seems the most up to date )
import pat3.utils as p3_u
import pat3.vehicles.fixed_wing.avl as p3_avl

def test0(avl_infile, avl_outdir, run_st=False, run_eig=False):
    P = p3_avl.make_pat_config(avl_infile, avl_outdir, run_st, run_eig)
    print(P)
    P.print_xml('/tmp/{}_avl.xml'.format(P.name))

    
def main():
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    avl_infile, out_dir = p3_u.pat_ressource('config/avl/cularis.avl'), '/tmp'
    test0(avl_infile, out_dir, run_st=False, run_eig=False)
    #test1()
    
if __name__ == "__main__":
    main()
