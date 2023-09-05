# Note for creating/building docker image for palom ([reference](https://micromamba-docker.readthedocs.io/en/latest/advanced_usage.html#advanced-usages))

1. Create reference env on micromamba's docker image

    ```bash
    # Run bash in micromamba docker image with bind volume for writing out env
    # lock file
    docker run -it --rm --platform linux/amd64 -v "$(pwd)":/data mambaorg/micromamba:1.4.9 bash
    ```

    ```bash
    # Manually install known deps in `palom` env 
    micromamba create -y -n palom python=3.10 "scikit-image<0.20" scikit-learn "zarr<2.15" tifffile imagecodecs matplotlib tqdm scipy dask numpy loguru=0.5.3 "ome-types>0.3" "pydantic<2" pint napari-lazy-openslide yamale fire termcolor dask-image -c conda-forge


    # Use `pip install --dry-run` to verify, would only expect to see `opencv`,
    # and `palom`
    micromamba activate palom
    python -m pip install --dry-run palom
    # output: Would install opencv-python-4.8.0.76 palom-2023.8.1


    # if the above checks out, export micromamba env as lock file
    micromamba env export --explicit > /data/docker-env.lock


    # pip install the rest of the packages, note: use `opencv-python-headless`
    # instead of `opencv-python`
    python -m pip install --no-deps palom==2023.8.1 opencv-python-headless==4.8.0.76


    # Test the environment
    python -c "import cv2; cv2.blur"
    ```

1. When building the docker image, specify `--platform linux/amd64`

    ```bash
    docker build --platform linux/amd64 --tag test-s3seg-large .
    ```