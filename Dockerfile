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

ARG CUDA=12.2.2
FROM nvidia/cuda:${CUDA}-cudnn8-runtime-ubuntu20.04
# FROM directive resets ARGS, so we specify again (the value is retained if
# previously set).
ARG CUDA

# Use bash to support string substitution.
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update --quiet \
    && apt-get install --no-install-recommends --yes --quiet \
       build-essential \
        cmake \
        cuda-command-line-tools-$(cut -f1,2 -d- <<< ${CUDA//./-}) \
        git \
        hmmer \
        kalign \
        tzdata \
        wget \
        aria2 \
        rsync \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get autoremove --yes \
    && apt-get clean

# Compile HHsuite from source.
RUN git clone --branch v3.3.0 --single-branch https://github.com/soedinglab/hh-suite.git /tmp/hh-suite \
    && mkdir /tmp/hh-suite/build \
    && pushd /tmp/hh-suite/build \
    && cmake -DCMAKE_INSTALL_PREFIX=/opt/hhsuite .. \
    && make --jobs 4 && make install \
    && ln -s /opt/hhsuite/bin/* /usr/bin \
    && popd \
    && rm -rf /tmp/hh-suite

# Install conda package manager.
RUN wget -q -P /tmp \
  https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash /tmp/Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda \
    && rm /tmp/Miniconda3-latest-Linux-x86_64.sh

# Install Conda packages.
ENV PATH="/opt/conda/bin:$PATH"
ENV LD_LIBRARY_PATH="/opt/conda/lib:$LD_LIBRARY_PATH"
ENV CONDA_PLUGINS_AUTO_ACCEPT_TOS=true
RUN echo $PATH \
    && conda install --quiet --yes conda==24.11.1 pip python=3.9.18 \
    && conda install --quiet --yes --channel nvidia cuda=${CUDA_VERSION} \
    && conda install --quiet --yes --channel conda-forge openmm=8.0.0 pdbfixer \
    && conda install --quiet --yes --channel bioconda nextflow==24.10.4 \
    && conda clean --all --force-pkgs-dirs --yes

COPY . /app/alphafold
RUN wget -q -P /app/alphafold/alphafold/common/ \
  https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt

RUN wget -q https://mmseqs.com/latest/mmseqs-linux-gpu.tar.gz \
    && tar xvfz mmseqs-linux-gpu.tar.gz -C /app/alphafold/ \
    && rm mmseqs-linux-gpu.tar.gz

ENV PATH="/app/alphafold/mmseqs/bin/:$PATH"

# Install pip packages.
RUN pip3 install --upgrade pip --no-cache-dir \
    && pip3 install -r /app/alphafold/docker/requirements.txt --no-cache-dir

# Add SETUID bit to the ldconfig binary so that non-root users can run it.
RUN chmod u+s /sbin/ldconfig.real

# Currently needed to avoid undefined_symbol error.
RUN ln -sf /usr/lib/x86_64-linux-gnu/libffi.so.7 /opt/conda/lib/libffi.so.7

ENV AF_CACHE="/app/alphafold"
CMD ["/bin/bash"]
