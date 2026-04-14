#!/usr/bin/env python
import json
import os
import pathlib
import shutil
import hashlib
import datetime

from absl import app
from absl import flags

from alphafold3.data.pipeline import _get_protein_templates
from alphafold3.data import msa_config
from alphafold3.common import folding_input

flags.DEFINE_list(
    'fasta_paths', None, 'Paths to FASTA files, each containing a prediction '
    'target that will be folded one after another. If a FASTA file contains '
    'multiple sequences, then it will be folded as a multimer. Paths should be '
    'separated by commas. All FASTA paths must have a unique basename as the '
    'basename is used to name the output directories for each prediction.')
flags.DEFINE_string('json_cache', None, 'Path to a cache directory for monomer .pkl features')
flags.DEFINE_string('output_dir', None, 'Path to a directory that will '
                    'store the results.')
flags.DEFINE_boolean('templates', True, 'Enable template search in monomer and multimer pipeline')

DB_DIR = flags.DEFINE_string(
    'db_dir',
    'alphafold3_data/',
    'Path to the directory containing the databases.',
)
_SEQRES_DATABASE_PATH = flags.DEFINE_string(
    'seqres_database_path',
    '${DB_DIR}/pdb_seqres_2022_09_28.fasta',
    'PDB sequence database path, used for template search.',
)
_HMMSEARCH_BINARY_PATH = flags.DEFINE_string(
    'hmmsearch_binary_path',
    shutil.which('hmmsearch'),
    'Path to the Hmmsearch binary.',
)
_HMMBUILD_BINARY_PATH = flags.DEFINE_string(
    'hmmbuild_binary_path',
    shutil.which('hmmbuild'),
    'Path to the Hmmbuild binary.',
)
_MAX_TEMPLATE_DATE = flags.DEFINE_string(
    'max_template_date',
    '2021-09-30',  # By default, use the date from the AlphaFold 3 paper.
    'Maximum template release date to consider. Format: YYYY-MM-DD. All '
    'templates released after this date will be ignored.',
)

FLAGS = flags.FLAGS

def SimpleFastaParser(handle):
    # taken from Bio.SeqIO.FastaIO
    # Skip any text before the first record (e.g. blank lines, comments)
    for line in handle:
        if line[0] == ">":
            title = line[1:].rstrip()
            break
    else:
        # no break encountered - probably an empty file
        return

    # Main logic
    # Note, remove trailing whitespace, and any internal spaces
    # (and any embedded \r which are possible in mangled files
    # when not opened in universal read lines mode)
    lines = []
    for line in handle:
        if line[0] == ">":
            yield title, "".join(lines).replace(" ", "").replace("\r", "")
            lines = []
            title = line[1:].rstrip()
            continue
        lines.append(line.rstrip())

    yield title, "".join(lines).replace(" ", "").replace("\r", "")


def main(_):

  _templates_config = msa_config.TemplatesConfig(
      template_tool_config=msa_config.TemplateToolConfig(
          database_path=FLAGS.seqres_database_path,
          chain_poly_type='polypeptide(L)',
          hmmsearch_config=msa_config.HmmsearchConfig(
              hmmsearch_binary_path=FLAGS.hmmsearch_binary_path,
              hmmbuild_binary_path=FLAGS.hmmbuild_binary_path,
              filter_f1=0.1,
              filter_f2=0.1,
              filter_f3=0.1,
              e_value=100,
              inc_e=100,
              dom_e=100,
              incdom_e=100,
              alphabet='amino',
          ),
      ),
      filter_config=msa_config.TemplateFilterConfig(
          max_subsequence_ratio=0.95,
          min_align_ratio=0.1,
          min_hit_length=10,
          deduplicate_sequences=True,
          max_hits=4,
          max_template_date=datetime.date.fromisoformat(FLAGS.max_template_date),
      ),
  )
  fasta_names = [pathlib.Path(p).stem for p in FLAGS.fasta_paths]

  for i, fasta_path in enumerate(FLAGS.fasta_paths):
    fasta_name = fasta_names[i]
    output_dir = os.path.join(FLAGS.output_dir, fasta_name)
    
    paired_msa_path = os.path.realpath(os.path.join(output_dir, "paired_hits.a3m"))
    unpaired_msa_path = os.path.realpath(os.path.join(output_dir, "unpaired_hits.a3m"))

    target_sequences_dict = []
    with open(fasta_path, 'r') as f:
        target_sequence = [record for record in SimpleFastaParser(f)][0][1]
    #target_sequence = [record.seq for record in SeqIO.parse(fasta_path, "fasta")][0]
    sequence_hash = hashlib.md5(str(target_sequence).encode()).hexdigest()
    out_json = os.path.join(FLAGS.json_cache, f"{sequence_hash}.json")

    if FLAGS.templates:
      template_hits = _get_protein_templates(
        sequence = target_sequence,
        input_msa_a3m = "".join(open(paired_msa_path, "r").readlines()),
        run_template_search = True,
        templates_config = _templates_config,
        pdb_database_path = f"{FLAGS.db_dir}/mmcif_files",
      )

      templates = [
            folding_input.Template(
                mmcif=struc.to_mmcif(),
                query_to_template_map=hit.query_to_hit_mapping,
            )
            for hit, struc in template_hits.get_hits_with_structures()
      ]
      ser_templates = [
            {
                'mmcif': template.mmcif,
                'queryIndices': list(template.query_to_template_map.keys()),
                'templateIndices': (
                    list(template.query_to_template_map.values()) or None
                ),
            }
            for template in templates
      ]
    else:
      ser_templates = []

    chain_dict = {"protein": {"id": "A",
                              "sequence": str(target_sequence),
                              "unpairedMsaPath": str(unpaired_msa_path),
                              "pairedMsaPath": str(paired_msa_path),
                              "templates": ser_templates,
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


if __name__ == '__main__':
  app.run(main)
