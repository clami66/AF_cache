# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Full AlphaFold protein structure prediction script."""
import enum
import json
import os
import pathlib
import shutil

from absl import app
from absl import flags
from absl import logging

flags.DEFINE_list(
    'fasta_paths', None, 'Paths to FASTA files, each containing a prediction '
    'target that will be folded one after another. If a FASTA file contains '
    'multiple sequences, then it will be folded as a multimer. Paths should be '
    'separated by commas. All FASTA paths must have a unique basename as the '
    'basename is used to name the output directories for each prediction.')
flags.DEFINE_string('json_cache', None, 'Path to a cache directory for monomer .pkl features')
flags.DEFINE_string('output_dir', None, 'Path to a directory that will '
                    'store the results.')
FLAGS = flags.FLAGS

def parse_af3():
  from Bio import SeqIO
  import hashlib
  fasta_names = [pathlib.Path(p).stem for p in FLAGS.fasta_paths]

  for i, fasta_path in enumerate(FLAGS.fasta_paths):
    fasta_name = fasta_names[i]
    output_dir = os.path.join(FLAGS.output_dir, fasta_name)
    
    paired_msa_path = os.path.realpath(os.path.join(output_dir, "paired_hits.a3m"))
    unpaired_msa_path = os.path.realpath(os.path.join(output_dir, "unpaired_hits.a3m"))

    target_sequences_dict = []
    target_sequence = [record.seq for record in SeqIO.parse(fasta_path, "fasta")][0]
    sequence_hash = hashlib.md5(target_sequence.encode()).hexdigest()
    out_json = os.path.join(FLAGS.json_cache, f"{sequence_hash}.json")
    chain_dict = {"protein": {"id": "A",
                              "sequence": str(target_sequence),
                              "unpairedMsaPath": str(unpaired_msa_path),
                              "pairedMsaPath": str(paired_msa_path),
                              "templates": None,
                              }}
    target_sequences_dict.append(chain_dict)
    json_data = {"name": "parse",
                 "modelSeeds": [1],
                 "sequences": target_sequences_dict,
                 "dialect": "alphafold3",
                 "version": 2,
                }
    with open(out_json, "w") as out:
      json.dump(json_data, out, indent=4, sort_keys=True)
  return out_json


def main(argv):
  parse_af3()


if __name__ == '__main__':
  app.run(main)
