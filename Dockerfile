FROM python:3.6
sudo apt update
sudo apt install libgl1-mesa-glx
RUN pip install scikit-learn scikit-image matplotlib tifffile opencv-python joblib
COPY S3segmenter.py ./app/S3segmenter.py
COPY rowit.py ./app/rowit.py
