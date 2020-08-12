FROM python:3.6

RUN pip install scikit-learn scikit-image==0.17.2 matplotlib tifffile==2020.7.24 opencv-python==4.4.0.40 joblib
sudo apt update
sudo apt install libgl1-mesa-glx
COPY S3segmenter.py ./app/S3segmenter.py
COPY rowit.py ./app/rowit.py
