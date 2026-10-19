#! /usr/bin/env python

from cupy import *
import sys

high_res = 3.2

par = Param(open("cns_xtal.par"))
pdbin = "minim_corrected.pdb"


#fix_pdb(pdbin, [["OH2 HOH Z","OH2 WAT Z"]])
par.coordinate_infile = pdbin
#par.ligand_topology_infile = "/home/legrand/lib/toppar/ligand.top"
#par.ligand_parameter_infile = "/home/legrand/lib/toppar/ligand.par"
par.ligand_topology_infile = ""
par.ligand_parameter_infile = ""
par.disulphide_dist=3.6
#run_cns("generate_easy.inp", par, "generate_easy.out")


par.structure_infile = "generate_easy.mtf"
par.parameter_infile_1 = "CNS_TOPPAR:protein_rep.param"
par.parameter_infile_2 = par.ligand_parameter_infile
par.parameter_infile_3 = "CNS_TOPPAR:water_rep.param"

par.high_res = 3.2
par.rigid_nstep = 10
par.num_cycles = 2
par.coordinate_infile = "generate_easy.pdb"
par.coordinate_outfile = "rigid_1.pdb"
par.r_group_1 = "(segid A )"
par.r_group_2 = "(segid B )"
par.r_group_3 = "(segid C )"
par.r_group_4 = "(segid D )"
par.bscale = "aniso"
#run_cns("rigid.inp", par, "rigid_1.out")

par.ncs_infile = "ncs.def"
par.high_res = high_res
par.coordinate_infile = "rigid_1.pdb"
par.coordinate_outfile = "minimize_y.pdb"
par.minimize_nstep=25;
par.num_cycles=3;
#run_cns("minimize.inp", par, "minimize_1.out")


n = 4
if 4 == 4:
    
    par.coordinate_infile = "minimize_y.pdb"
    par.output_root = "refine_%d" % (n+1)
    par.anneal = "true"
    par.atom_rigid = "(segid N or segid M or segid S)"
    par.md_scheme = "slowcool"
    par.temperature = 1200
    par.cool_rate = 20
    par.num_cycles=3
    par.minimize_nstep = 60
    par.torsion_maxtree = 20
    par.torsion_maxchain = 3000
    par.bfactor_nstep = 10
    run_cns("refine.inp", par, "refine_%d.out" % (n+1))

sys.exit()

par.coordinate_infile = par.coordinate_outfile
par.output_root = "refine_1"
par.anneal = "false"
par.md_scheme = "slowcool"
par.temperature = 3000
par.cool_rate = 50
par.num_cycles=1
par.minimize_nstep = 30
par.bfactor_nstep = 20
run_cns("refine.inp", par, "refine_1.out")

for n in range(1,4):
    
    par.coordinate_infile = "refine_%d.pdb" % n
    par.coordinate_outfile = "water_pick_%d.pdb" % n
    par.structure_outfile = "water_pick_%d.mtf" % n
    par.wat_resname="WAT"
    par.wat_segid="W"
    par.peak=3.0
    par.minimize_nstep = 20
    par.bfactor_nstep = 20
    run_cns("water_pick.inp", par, "water_pick_%d.out" % n)
    
    par.structure_infile  = par.structure_outfile
    par.coordinate_infile = par.coordinate_outfile
    par.output_root = "refine_%d" % (n+1)
    par.anneal = "true"
    par.atom_rigid = "(segid N or segid M or segid S)"
    par.md_scheme = "slowcool"
    par.temperature = 3000
    par.cool_rate = 50
    par.num_cycles=1
    par.minimize_nstep = 30
    par.torsion_maxtree = 10
    par.bfactor_nstep = 20
    run_cns("refine.inp", par, "refine_%d.out" % (n+1))


par.structure_infile  = "generate_easy.mtf"
par.coordinate_infile = "refine_3.pdb"
par.output_root = "anneal_4"
par.atom_fixed = "(segid N)"
par.temperature = 4000
par.cool_rate = 50
par.num_trials = 4
par.final_dyn = 100
par.final_min = 50
run_cns("anneal.inp", par, "anneal.out")

par.structure_infile  = "generate_easy.mtf"
par.coordinate_infile = par.coordinate_outfile 
par.output_root = "anneal_4"
par.atom_fixed = "(segid N or segid M)"
par.temperature = 5000
par.cool_rate = 50
par.num_trials = 4
par.final_dyn = 100
par.final_min = 50
run_cns("xtal_refine/anneal.inp", par, "anneal.out")
