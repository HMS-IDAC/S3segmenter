FROM mambaorg/micromamba:1.4.9

COPY --chown=$MAMBA_USER:$MAMBA_USER docker-env.lock /tmp/docker-env.lock
RUN micromamba install --name base --yes --file /tmp/docker-env.lock \
    && micromamba clean --trash -aflp --yes

# pip install packages that are not available/problematic on conda-forge
RUN /opt/conda/bin/python -m pip install \
    --no-deps \
    opencv-python-headless==4.8.0.76 \
    palom==2023.8.1 \
    && /opt/conda/bin/python -m pip cache purge

# copy scripts to /app/
COPY / /app/

# add conda path to PATH to allow entrypoint overwrite
ENV PATH="${PATH}:/opt/conda/bin"