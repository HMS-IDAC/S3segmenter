FROM python:3.10

RUN pip install scikit-learn scikit-image opencv-python tifffile imagecodecs dask-image zarr ome_types

COPY S3segmenter.py ./app/S3segmenter.py
COPY save_tifffile_pyramid.py ./app/save_tifffile_pyramid.py
COPY watershed.py ./app/watershed.py
