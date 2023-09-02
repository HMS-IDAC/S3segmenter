FROM python:3.10

RUN python -m pip install \
    palom>=2023.8.1 \ 
    dask[dataframe] \
    dask-image \
    ome_types \
    pydantic<2 \
    && python -m pip cache purge
RUN python -m pip uninstall opencv-python -y
RUN python -m pip install opencv-python-headless \
    && python -m pip cache purge

COPY S3segmenter.py ./app/S3segmenter.py
COPY save_tifffile_pyramid.py ./app/save_tifffile_pyramid.py
COPY watershed.py ./app/watershed.py
COPY s3seg_qc.py ./app/s3seg_qc.py
COPY s3seg_util.py ./app/s3seg_util.py
COPY ignored_args.py ./app/ignored_args.py

