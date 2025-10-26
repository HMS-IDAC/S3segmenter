FROM python:3.13.9

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y libgl1 \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache numpy scikit-learn scikit-image matplotlib tifffile opencv-python ome_types imagecodecs

COPY S3segmenter.py ./app/S3segmenter.py
COPY save_tifffile_pyramid.py ./app/save_tifffile_pyramid.py
COPY rowit.py ./app/rowit.py

