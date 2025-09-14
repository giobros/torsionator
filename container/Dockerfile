#FROM nvidia/cuda:11.8.0-cudnn8-runtime-ubuntu20.04
FROM nvidia/cuda:11.8.0-cudnn8-devel-ubuntu20.04
ENV LD_LIBRARY_PATH="/usr/local/cuda-11.8/lib64":$LD_LIBRARY_PATH
ENV DEBIAN_FRONTEND=noninteractive
ENV OBI_ROOT="/app/OBIWAN"
ENV OBIWAN_MODEL_PATH="${OBI_ROOT}/results/models"
ENV PATH="/opt/conda/bin:$PATH"

# ------------------------------------
# 1. Install system dependencies
# ------------------------------------
RUN apt-get update && \
    apt-get install -y wget curl git build-essential && \
    rm -rf /var/lib/apt/lists/*

# ------------------------------------
# 2. Install Miniconda + Mamba
# ------------------------------------
# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    /opt/conda/bin/conda config --system --remove channels defaults && \
    /opt/conda/bin/conda config --system --add channels conda-forge && \
    /opt/conda/bin/conda install -y -n base mamba
# ------------------------------------
# 3. Create env using Mamba
# ------------------------------------
RUN /opt/conda/bin/mamba create -y -n obiwan_env -c conda-forge python=3.9 ambertools=22

# ------------------------------------
# 4. Install Python packages
# ------------------------------------
RUN /opt/conda/bin/mamba run -n obiwan_env pip install --no-cache-dir \
    tensorflow==2.14 \
    numpy==1.23.5 \
    Pillow \
    matplotlib \
    scipy \
    pandas==2.0.3 \
    tqdm \
    rdkit==2023.03.3 \
    joblib \
    networkx 

RUN /opt/conda/bin/mamba run -n obiwan_env pip install --no-cache-dir torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118 
RUN /opt/conda/bin/mamba run -n obiwan_env pip install  mace-torch
# ------------------------------------
# 5. Clone and install OBIWAN
# ------------------------------------
RUN git clone https://github.com/virtualmartire/OBIWAN.git ${OBI_ROOT}
WORKDIR ${OBI_ROOT}

# ------------------------------------
# 6. Copy custom ASE calculator
# ------------------------------------
COPY obi.py /opt/conda/envs/obiwan_env/lib/python3.9/site-packages/ase/calculators/obi.py
#COPY torsionator /opt/torsionator 
#RUN /opt/conda/envs/obiwan_env/bin/pip install -e /opt/torsionator
# ------------------------------------
# 7. Use mamba env by default
# ------------------------------------
ENV PATH="/opt/conda/envs/obiwan_env/bin:$PATH"
ENV CONDA_DEFAULT_ENV=obiwan_env

# ------------------------------------
# 8. Final working directory
# ------------------------------------
WORKDIR /workspace
