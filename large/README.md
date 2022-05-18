TODO
- write QC files
- parse and write pixel size to resulting masks

---

conda env for development and testing

```bash
conda create -n s3seg-dev numpy pandas scikit-learn scikit-image imagecodecs tifffile zarr dask-image -c conda-forge
conda activate s3seg-dev
python -m pip install ome_types opencv-python
```

---

input file

```
Y:\SORGER\DATA\COMPUTATION\YU-AN\TB_B5
    TB_B5_1120_RCPNL-pmap.tif
```


output files

```
Y:\SORGER\DATA\COMPUTATION\YU-AN\TB_B5
│   TB_B5_1120_RCPNL-pmap.tif
│
└───segmentation
    └───TB_B5_1120_RCPNL
            cell.ome.tif
            cyto.ome.tif
            nuclei.ome.tif
```


run script 

```bash 
python S3segmenter.py \
    --imagePath "TB_B5_1120_RCPNL" \
    --stackProbPath "Y:\sorger\data\computation\Yu-An\TB_B5\TB_B5_1120_RCPNL-pmap.tif" \
    --outputPath "Y:\sorger\data\computation\Yu-An\TB_B5\segmentation" \
    --area-max 5000 --expand-size 3 --maxima-footprint-size 5 --mean-intensity-min 80 
```


console log

```
2022-05-12 04:25:54 | INFO     | Reading Y:\sorger\data\computation\Yu-An\TB_B5\TB_B5_1120_RCPNL-pmap.tif (watershed.py:388)
2022-05-12 04:26:19 | INFO     | Probability map shape: (2, 33045, 42992) (watershed.py:390)
2022-05-12 04:26:22 | INFO     | Run initiating (watershed.py:102)
2022-05-12 04:26:22 | INFO     | Configuration:

{'area_max': 5000,
 'area_min': 10,
 'footprint_size': 5,
 'h': 0.01,
 'mean_intensity_min': 80.0,
 'overlap_size': 128,
 'sigma': 1}

 (watershed.py:103)
2022-05-12 04:26:22 | INFO     | Start configuring `label` (watershed.py:196)
2022-05-12 04:26:30 | INFO     | End configuring `label` (watershed.py:198)
2022-05-12 04:26:30 | INFO     | Start configuring `label` (watershed.py:109)
2022-05-12 04:27:02 | INFO     | End configuring `label` (watershed.py:111)
[########################################] | 100% Completed |  9min 52.0s
2022-05-12 04:36:57 | INFO     | Writing to Y:\sorger\data\computation\Yu-An\TB_B5\segmentation\TB_B5_1120_RCPNL\nuclei.ome.tif (watershed.py:127)
2022-05-12 04:38:28 | INFO     | Expanding 3 pixels (watershed.py:411)
2022-05-12 04:38:28 | INFO     | Writing to Y:\sorger\data\computation\Yu-An\TB_B5\segmentation\TB_B5_1120_RCPNL\cell.ome.tif (watershed.py:127)
[########################################] | 100% Completed |  1min 13.2s
[########################################] | 100% Completed | 26.5s
2022-05-12 04:41:40 | INFO     | Generating difference mask (watershed.py:422)
2022-05-12 04:41:40 | INFO     | Writing to Y:\sorger\data\computation\Yu-An\TB_B5\segmentation\TB_B5_1120_RCPNL\cyto.ome.tif (watershed.py:127)
[########################################] | 100% Completed |  1min  3.1s
[########################################] | 100% Completed |  7.7s
2022-05-12 04:44:24 | INFO     | Done (watershed.py:433)
```
