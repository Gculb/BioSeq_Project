FROM mambaorg/micromamba:2.0.5

COPY --chown=$MAMBA_USER:$MAMBA_USER containers/environment.yml /tmp/environment.yml
RUN micromamba create --yes --name seqcheckflow --file /tmp/environment.yml && \
    micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER seqcheckflow /workspace/seqcheckflow

ENV PATH=/opt/conda/envs/seqcheckflow/bin:$PATH
ENV PYTHONPATH=/workspace
WORKDIR /workspace
