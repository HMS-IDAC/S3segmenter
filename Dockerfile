FROM python:3.10.2

RUN pip install numpy==1.22.3 scikit-learn scikit-image matplotlib tifffile opencv-python ome_types

COPY S3segmenter.py ./app/S3segmenter.py
COPY save_tifffile_pyramid.py ./app/save_tifffile_pyramid.py
COPY rowit.py ./app/rowit.py

