FROM python:3.6

RUN pip install scikit-learn scikit-image matplotlib tifffile opencv-python joblib
sudo apt update
sudo apt install libgl1-mesa-glx
COPY S3segmenter.py ./app/S3segmenter.py
COPY rowit.py ./app/rowit.py
