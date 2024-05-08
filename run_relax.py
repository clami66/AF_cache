import json
import os
import pathlib
import pickle
import random
import sys
import time
import argparse
from typing import Dict

from absl import app
from absl import flags
from absl import logging
from alphafold.common import protein
from alphafold.common import residue_constants
from alphafold.data import pipeline
from alphafold.data import templates
from alphafold.model import data
from alphafold.model import config
from alphafold.model import model
from alphafold.relax import relax
import numpy as np

logging.set_verbosity(logging.INFO)
logging.use_absl_handler()

parser = argparse.ArgumentParser(description='Run Amber Relax (AlphaFold2 settings) on any structure')
parser.add_argument('model_name', help = 'PDB File to Relax')
parser.add_argument('--output_dir', type = str, help = 'Output path for the relaxed model')
parser.add_argument('--use_gpu', action = 'store_true', help = 'Whether relaxing should be done with GPU (default is CPU)')
# An option to tell where the AF2 code is. 
# Also change the **default** to wherever your cloned AF2 is.
parser.add_argument("--af2_dir", default="**your/path/to/alphafold-2.3.1/**", help="AlphaFold code directory")

args = parser.parse_args()
# This line below should put the AF code in your environment path
sys.path.append(args.af2_dir)
model  = args.model_name

def relax_with_amber(model_name,output_dir):
        MAX_TEMPLATE_HITS = 20
        RELAX_MAX_ITERATIONS = 0
        RELAX_ENERGY_TOLERANCE = 2.39
        RELAX_STIFFNESS = 10.0
        RELAX_EXCLUDE_RESIDUES = []
        RELAX_MAX_OUTER_ITERATIONS = 3
        amber_relax = relax.AmberRelaxation(max_iterations=RELAX_MAX_ITERATIONS, tolerance=RELAX_ENERGY_TOLERANCE, stiffness=RELAX_STIFFNESS, exclude_residues=RELAX_EXCLUDE_RESIDUES, max_outer_iterations=RELAX_MAX_OUTER_ITERATIONS, use_gpu=args.use_gpu)
        unrelaxed_protein = model_name

        logging.info("Relaxing model %s", unrelaxed_protein)
        with open(str(model_name)) as f:
                test_prot = protein.from_pdb_string(f.read())
                pdb_min, _, _ = amber_relax.process(prot=test_prot)
                print(pdb_min)
                with open(str(output_dir)+'amber_r_'+str(model_name).split('/')[-1],'w+') as rel_f:
                        rel_f.write(str(pdb_min))

relax_with_amber(model,args.output_dir)
