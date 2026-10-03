FROM mambaorg/micromamba:2.0.5

COPY --chown=$MAMBA_USER:$MAMBA_USER containers/environment.yml /tmp/environment.yml
RUN micromamba create --yes --name bioseq --file /tmp/environment.yml && \
    micromamba clean --all --yes

ENV PATH=/opt/conda/envs/bioseq/bin:$PATH
ENV PYTHONPATH=/workspace
WORKDIR /workspace
