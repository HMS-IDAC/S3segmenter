FROM python:3.6

RUN pip install scikit-learn scikit-image matplotlib tifffile opencv-python

COPY S3segmenter.py ./app/S3segmenter.py
