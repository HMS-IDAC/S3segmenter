FROM python:3.8.0

RUN pip install numpy==1.22.3 scikit-learn scikit-image==0.16.2 matplotlib tifffile==2021.6.6 opencv-python==4.3.0.36 ome_types==0.4.5 imagecodecs

COPY S3segmenter.py ./app/S3segmenter.py
COPY save_tifffile_pyramid.py ./app/save_tifffile_pyramid.py
COPY rowit.py ./app/rowit.py

