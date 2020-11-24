#! /usr/bin/env python
import logging, os, math, numpy as np
import pdb
#
# Note to self:
# /mnt/mint17/home/poine/work/python-aerospace-toolbox/config/avl/cularis.avl
# 
# /mnt/mint17/home/poine/work/ngfw/trunk/simulations/py_avl/config/avl/cularis.avl
# ( seems the most up to date )
#
import pat3.utils as p3_u
import pat3.vehicles.fixed_wing.avl as p3_avl

def test0(avl_infile, avl_outdir, force_run_st=False, force_run_eig=False):
    P = p3_avl.make_pat_config(avl_infile, avl_outdir, force_run_st, force_run_eig)
    print(P)
    P.save_xml('/tmp/{}_avl.xml'.format(P.name))

    
def main():
    logging.basicConfig(level=logging.INFO); np.set_printoptions(linewidth=500)
    for _ac in ['cularis']:#, 'funjet']:
        avl_infile, out_dir = p3_u.pat_ressource('config/avl/{}.avl'.format(_ac)), '/tmp'
        test0(avl_infile, out_dir, force_run_st=False, force_run_eig=False)


if __name__ == "__main__":
    main()
