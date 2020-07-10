FROM python:3.6

RUN pip install scikit-learn scikit-image==0.17.2 matplotlib tifffile==2020.7.4 opencv-python

COPY S3segmenter.py ./app/S3segmenter.py
