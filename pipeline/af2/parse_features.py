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

logging.set_verbosity(logging.INFO)
MAX_TEMPLATE_HITS = 20

flags.DEFINE_list(
    'fasta_paths', None, 'Paths to FASTA files, each containing a prediction '
    'target that will be folded one after another. If a FASTA file contains '
    'multiple sequences, then it will be folded as a multimer. Paths should be '
    'separated by commas. All FASTA paths must have a unique basename as the '
    'basename is used to name the output directories for each prediction.')

flags.DEFINE_string('output_dir', None, 'Path to a directory that will '
                    'store the results.')
flags.DEFINE_string('jackhmmer_binary_path', shutil.which('jackhmmer'),
                    'Path to the JackHMMER executable.')
flags.DEFINE_string('hhblits_binary_path', shutil.which('hhblits'),
                    'Path to the HHblits executable.')
flags.DEFINE_string('hhsearch_binary_path', shutil.which('hhsearch'),
                    'Path to the HHsearch executable.')
flags.DEFINE_string('hmmsearch_binary_path', shutil.which('hmmsearch'),
                    'Path to the hmmsearch executable.')
flags.DEFINE_string('hmmbuild_binary_path', shutil.which('hmmbuild'),
                    'Path to the hmmbuild executable.')
flags.DEFINE_string('kalign_binary_path', shutil.which('kalign'),
                    'Path to the Kalign executable.')
flags.DEFINE_string('mmseqs2_binary_path', shutil.which('mmseqs'),
                    'Path to the MMseqs2 executable.')
flags.DEFINE_string('mmseqs2_uniref_database_path', None, 'Path to the Uniref30 '
                    'database for use by MMseqs2.')
flags.DEFINE_string('mmseqs2_env_database_path', None, 'Path to the environmental '
                    'database for use by MMseqs2.')
flags.DEFINE_string('uniref90_database_path', None, 'Path to the Uniref90 '
                    'database for use by JackHMMER.')
flags.DEFINE_string('mgnify_database_path', None, 'Path to the MGnify '
                    'database for use by JackHMMER.')
flags.DEFINE_string('bfd_database_path', None, 'Path to the BFD '
                    'database for use by HHblits.')
flags.DEFINE_string('small_bfd_database_path', None, 'Path to the small '
                    'version of BFD used with the "reduced_dbs" preset.')
flags.DEFINE_string('uniref30_database_path', None, 'Path to the UniRef30 '
                    'database for use by HHblits.')
flags.DEFINE_string('uniprot_database_path', None, 'Path to the Uniprot '
                    'database for use by JackHMMer.')
flags.DEFINE_string('pdb70_database_path', None, 'Path to the PDB70 '
                    'database for use by HHsearch.')
flags.DEFINE_string('pdb_seqres_database_path', None, 'Path to the PDB '
                    'seqres database for use by hmmsearch.')
flags.DEFINE_string('template_mmcif_dir', None, 'Path to a directory with '
                    'template mmCIF structures, each named <pdb_id>.cif')
flags.DEFINE_string('max_template_date', None, 'Maximum template release date '
                    'to consider. Important if folding historical test sets.')
flags.DEFINE_string('obsolete_pdbs_path', None, 'Path to file containing a '
                    'mapping from obsolete PDB IDs to the PDB IDs of their '
                    'replacements.')
flags.DEFINE_string('pickle_cache', None, 'Path to a cache directory for monomer .pkl features')
flags.DEFINE_enum('db_preset', 'full_dbs',
                  ['full_dbs', 'reduced_dbs'],
                  'Choose preset MSA database configuration - '
                  'smaller genetic database config (reduced_dbs) or '
                  'full genetic database config  (full_dbs)')
flags.DEFINE_enum('model_preset', 'monomer',
                  ['monomer', 'monomer_casp14', 'monomer_ptm', 'multimer', 'multimer_v1', 'multimer_v2', 'multimer_v3'],
                  'Choose preset model configuration - the monomer model, '
                  'the monomer model with extra ensembling, monomer model with '
                  'pTM head, or multimer model')
flags.DEFINE_boolean('use_precomputed_msas', True, 'Whether to read MSAs that '
                     'have been written to disk instead of running the MSA '
                     'tools. The MSA files are looked up in the output '
                     'directory, so it must stay the same between multiple '
                     'runs that are to reuse the MSAs. WARNING: This will not '
                     'check if the sequence, database or configuration have '
                     'changed.')
flags.DEFINE_integer('uniprot_max_hits', 50000, 'Max hits in uniprot MSA')
flags.DEFINE_integer('mgnify_max_hits', 500, 'Max hits in uniprot MSA')
flags.DEFINE_integer('uniref_max_hits', 10000, 'Max hits in uniprot MSA')
flags.DEFINE_integer('bfd_max_hits', 10000, 'Max hits in uniprot MSA')
flags.DEFINE_boolean('separate_homomer_msas', False, 'Whether to force separate processing of homomer MSAs')
flags.DEFINE_integer('redundancy_reduce_templates', 100, 'Percentage of redundancy reduction for template hits')
flags.DEFINE_list('pad_to_size', [None, None], 'Pad input features to a given seq. length x MSA depth. '
                     'This is useful when processing multiple sequences at once to avoid re-compiling the models')

FLAGS = flags.FLAGS

def _check_flag(flag_name: str,
                other_flag_name: str,
                should_be_set: bool):
  if should_be_set != bool(FLAGS[flag_name].value):
    verb = 'be' if should_be_set else 'not be'
    raise ValueError(f'{flag_name} must {verb} set when running with '
                     f'"--{other_flag_name}={FLAGS[other_flag_name].value}".')

def main(argv):
  # delay imports
  from alphafold.data import pipeline
  from alphafold.data import pipeline_multimer
  from alphafold.data import templates
  from alphafold.data.tools import hhsearch
  from alphafold.data.tools import hmmsearch
  from alphafold.model import config

  for tool_name in (
      'jackhmmer', 'hhblits', 'hhsearch', 'hmmsearch', 'hmmbuild', 'kalign'):
    if not FLAGS[f'{tool_name}_binary_path'].value:
      raise ValueError(f'Could not find path to the "{tool_name}" binary. Make '
                      'sure it is installed on your system.')

  use_small_bfd = FLAGS.db_preset == 'reduced_dbs'
  _check_flag('small_bfd_database_path', 'db_preset',
              should_be_set=use_small_bfd)
  _check_flag('bfd_database_path', 'db_preset',
              should_be_set=not use_small_bfd)
  _check_flag('uniref30_database_path', 'db_preset',
              should_be_set=not use_small_bfd)

  run_multimer_system = 'multimer' in FLAGS.model_preset
  _check_flag('pdb70_database_path', 'model_preset',
              should_be_set=not run_multimer_system)
  _check_flag('pdb_seqres_database_path', 'model_preset',
              should_be_set=run_multimer_system)
  _check_flag('uniprot_database_path', 'model_preset',
              should_be_set=run_multimer_system)

  # Check for duplicate FASTA file names.
  fasta_names = [pathlib.Path(p).stem for p in FLAGS.fasta_paths]
  if len(fasta_names) != len(set(fasta_names)):
    raise ValueError('All FASTA paths must have a unique basename.')

  if run_multimer_system:
    logging.info("redundancy_reduce FLAG: %d", FLAGS.redundancy_reduce_templates)
    template_searcher = hmmsearch.Hmmsearch(
        binary_path=FLAGS.hmmsearch_binary_path,
        hmmbuild_binary_path=FLAGS.hmmbuild_binary_path,
        database_path=FLAGS.pdb_seqres_database_path,
        redundancy_reduce=FLAGS.redundancy_reduce_templates)
    template_featurizer = templates.HmmsearchHitFeaturizer(
        mmcif_dir=FLAGS.template_mmcif_dir,
        max_template_date=FLAGS.max_template_date,
        max_hits=MAX_TEMPLATE_HITS,
        kalign_binary_path=FLAGS.kalign_binary_path,
        release_dates_path=None,
        obsolete_pdbs_path=FLAGS.obsolete_pdbs_path)
  else:
    template_searcher = hhsearch.HHSearch(
        binary_path=FLAGS.hhsearch_binary_path,
        databases=[FLAGS.pdb70_database_path])
    template_featurizer = templates.HhsearchHitFeaturizer(
        mmcif_dir=FLAGS.template_mmcif_dir,
        max_template_date=FLAGS.max_template_date,
        max_hits=MAX_TEMPLATE_HITS,
        kalign_binary_path=FLAGS.kalign_binary_path,
        release_dates_path=None,
        obsolete_pdbs_path=FLAGS.obsolete_pdbs_path)

  monomer_data_pipeline = pipeline.DataPipeline(
      jackhmmer_binary_path=FLAGS.jackhmmer_binary_path,
      hhblits_binary_path=FLAGS.hhblits_binary_path,
      mmseqs2_binary_path="mmseqs",
      uniref90_database_path=FLAGS.uniref90_database_path,
      mgnify_database_path=FLAGS.mgnify_database_path,
      mmseqs2_uniref_database_path=FLAGS.mmseqs2_uniref_database_path,
      mmseqs2_env_database_path=FLAGS.mmseqs2_env_database_path,
      bfd_database_path=FLAGS.bfd_database_path,
      uniref30_database_path=FLAGS.uniref30_database_path,
      small_bfd_database_path=FLAGS.small_bfd_database_path,
      template_searcher=template_searcher,
      template_featurizer=template_featurizer,
      use_small_bfd=use_small_bfd,
      use_precomputed_msas=FLAGS.use_precomputed_msas,
      mgnify_max_hits=FLAGS.mgnify_max_hits,
      uniref_max_hits=FLAGS.uniref_max_hits,
      bfd_max_hits=FLAGS.bfd_max_hits,
      alignments_only=True,
      no_uniref=True,
      no_mgnify=True)

  if run_multimer_system:
    data_pipeline = pipeline_multimer.DataPipeline(
        monomer_data_pipeline=monomer_data_pipeline,
        jackhmmer_binary_path=FLAGS.jackhmmer_binary_path,
        uniprot_database_path=FLAGS.uniprot_database_path,
        use_mmseqs2_align=True,
        use_precomputed_msas=FLAGS.use_precomputed_msas,
        max_uniprot_hits=FLAGS.uniprot_max_hits,
        separate_homomer_msas=FLAGS.separate_homomer_msas,
        feat_config=config.model_config("model_1_multimer_v3").eval.feat,
        pad_length=[int(FLAGS.pad_to_size[0]), int(FLAGS.pad_to_size[1])] if FLAGS.pad_to_size[0] else None,
        pickle_cache=FLAGS.pickle_cache)
  else:
    data_pipeline = monomer_data_pipeline

  # Predict structure for each of the sequences.
  for i, fasta_path in enumerate(FLAGS.fasta_paths):
    fasta_name = fasta_names[i]
    output_dir = os.path.join(FLAGS.output_dir, fasta_name)
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)
    msa_output_dir = os.path.join(output_dir, 'msas')
    if not os.path.exists(msa_output_dir):
      os.makedirs(msa_output_dir)
    logging.info("Parsing the features")
    # Get features.
    feature_dict = data_pipeline.process(
        input_fasta_path=fasta_path,
        msa_output_dir=msa_output_dir)


if __name__ == '__main__':
  flags.mark_flags_as_required([
      'fasta_paths',
      'output_dir',
      'uniref90_database_path',
      'mgnify_database_path',
      'template_mmcif_dir',
      'max_template_date',
      'obsolete_pdbs_path',
  ])

  app.run(main)
