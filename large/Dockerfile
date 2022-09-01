FROM python:3.10

RUN python -m pip install \
    palom \ 
    dask[dataframe] \
    opencv-python-headless \
    dask-image \
    ome_types

COPY S3segmenter.py ./app/S3segmenter.py
COPY save_tifffile_pyramid.py ./app/save_tifffile_pyramid.py
COPY watershed.py ./app/watershed.py
COPY s3seg_qc.py ./app/s3seg_qc.py
COPY s3seg_util.py ./app/s3seg_util.py

