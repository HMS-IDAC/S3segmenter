FROM python:3.6

RUN pip install scikit-learn scikit-image==0.16.2 matplotlib tifffile==0.15.1 opencv-python

COPY S3segmenter.py ./app/S3segmenter.py
