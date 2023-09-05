import pathlib
import numpy as np
import skimage.segmentation
import skimage.exposure
import dask.array as da
import zarr
import tifffile
import palom

import _version


def mask_to_bound(
    img_da,
    overlap_depth=128,
    out_dtype=np.uint8,
    positive_value=255
):
    return img_da.map_overlap(
        skimage.segmentation.find_boundaries,
        mode='outer',
        depth=overlap_depth,
        boundary='none'
    ).astype(out_dtype) * positive_value


def rescale_channel(
    img_da,
    min_percentile=0,
    max_percentile=100,
    out_dtype=np.uint8
):
    assert img_da.ndim == 2
    rescale_min = da.percentile(img_da.flatten(), min_percentile).compute()[0]
    rescale_max = da.percentile(img_da.flatten(), max_percentile).compute()[0]
    return img_da.map_blocks(
        skimage.exposure.rescale_intensity,
        in_range=(rescale_min, rescale_max),
        out_range=out_dtype,
        dtype=out_dtype
    ).astype(out_dtype)


def run_mcmicro(
    mask_path,
    out_path,
    pmap_path=None,
    img_path=None,
    img_channels=None,
    pixel_size=1
):
    assert pathlib.Path(out_path).parent.exists()
    
    out_channels = []
    channel_names = []

    zarr_mask = zarr.open( 
        tifffile.imread(mask_path, aszarr=True, series=0, level=0) 
    )
    da_mask = da.from_zarr(zarr_mask)
    da_outline = mask_to_bound(da_mask)
    out_channels.append(da_outline)
    channel_names.append('Mask inner outline')
    
    if pmap_path is not None:
        contour_probability_map = (
            palom.reader.OmePyramidReader(pmap_path).pyramid[0][1].compute()
        )
        da_contour_probability_map = da.from_array(
            contour_probability_map, chunks=2048
        )
        out_channels.append(da_contour_probability_map)
        channel_names.append('Contour probability map')
    
    if img_channels is None: img_channels = [0]
    img_channels = set(img_channels)

    if img_path is not None:
        reader = palom.reader.OmePyramidReader(img_path)
        da_img = reader.pyramid[0]

        for channel in img_channels:
            out_channels.append(
                rescale_channel(da_img[channel], 0.1, 99.9, np.uint8)
            )
            channel_names.append(f"Image channel {channel}")
    
    da_stack = da.array(out_channels)

    palom.pyramid.write_pyramid(
        [da_stack],
        out_path,
        channel_names=[channel_names],
        pixel_size=pixel_size,
        downscale_factor=2,
        compression='zlib',
        tile_size=1024,
        save_RAM=True,
        kwargs_tifffile=dict(software=f"s3segmenter-large v{_version.VERSION}")       
    )

    return 0
