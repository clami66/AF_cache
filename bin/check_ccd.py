#!/usr/bin/env python
import os
from alphafold3.common import resources
from alphafold3.common import safe_pickle
import subprocess

_CCD_SETS_CCD_PICKLE_FILE = resources.filename(
    resources.ROOT / 'constants/converters/chemical_component_sets.pickle'
)

print("Looking for CCD...")
if not os.path.exists(_CCD_SETS_CCD_PICKLE_FILE):
    print("Installing CCD...")
    subprocess.run(["build_data"])

print("Done.")

