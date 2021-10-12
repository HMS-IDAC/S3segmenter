FROM python:3.7.11

RUN pip install numpy==1.19 scikit-learn scikit-image==0.14.2 matplotlib tifffile==2021.6.6 opencv-python==4.3.0.36 ome_types

COPY S3segmenter.py ./app/S3segmenter.py
COPY save_tifffile_pyramid.py ./app/save_tifffile_pyramid.py
COPY rowit.py ./app/rowit.py

