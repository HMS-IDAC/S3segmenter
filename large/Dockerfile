FROM python:3.10

RUN python -m pip install \
    dask[dataframe] \
    scikit-learn \
    scikit-image \
    opencv-python-headless \
    tifffile \
    imagecodecs \
    dask-image \
    zarr \
    ome_types

COPY S3segmenter.py ./app/S3segmenter.py
COPY save_tifffile_pyramid.py ./app/save_tifffile_pyramid.py
COPY watershed.py ./app/watershed.py
