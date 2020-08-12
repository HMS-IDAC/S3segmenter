FROM ubuntu:19.10

RUN apt-get update \ 
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    python3.6 \
    python3-pip \
    libglu1 

RUN pip3 install -U \
	pip \
	scikit-image \
	matplotlib \
	tifffile \
	joblib \
	opencv-python
COPY S3segmenter.py ./app/S3segmenter.py
COPY rowit.py ./app/rowit.py

