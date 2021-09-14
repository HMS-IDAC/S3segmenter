import numpy as np
import tifffile
import skimage.transform

PHYSICAL_SIZE_UNIT = ['Ym', 'Zm', 'Em', 'Pm', 'Tm', 'Gm', 'Mm', 'km', 'hm', 'dam', 'm', 'dm', 'cm', 'mm', 'µm', 'nm', 'pm', 'fm', 'am', 'zm', 'ym', 'Å', 'thou', 'li', 'in', 'ft', 'yd', 'mi', 'ua', 'ly', 'pc', 'pt', 'pixel', 'reference frame']

def normalize_image_shape(img):
    assert img.ndim in (2, 3), (
        'image must be 2D (Y, X) or 3D (C, Y, X)'
    )
    if img.ndim == 2:
        img = img.reshape(1, *img.shape)
    assert np.argmin(img.shape) == 0, (
        '3D image must be in (C, Y, X) order'
    )
    return img

def save_pyramid(
    out_img, output_path,
    pixel_sizes=(1, 1),
    pixel_size_units=('µm', 'µm'),
    channel_names=None,
    software=None,
    is_mask=False
):
    assert '.ome.tif' in str(output_path)
    assert len(pixel_sizes) == len(pixel_size_units) == 2
    assert out_img.ndim in (2, 3), (
        'image must be either 2D (Y, X) or 3D (C, Y, X)'
    )
    
    img_shape_ori = out_img.shape
    out_img = normalize_image_shape(out_img)
    img_shape = out_img.shape

    size_x, size_y = np.array(pixel_sizes, dtype=float)
    unit_x, unit_y = pixel_size_units

    assert (unit_x in PHYSICAL_SIZE_UNIT) & (unit_y in PHYSICAL_SIZE_UNIT), (
        f'pixel_size_units must be a tuple of the followings: '
        f'{", ".join(PHYSICAL_SIZE_UNIT)}'
    )

    n_channels = img_shape[0]
    if channel_names == None:
        channel_names = [f'Channel {i}' for i in range(n_channels)]
    else:
        if type(channel_names) == str:
            channel_names = [channel_names]
        n_channel_names = len(channel_names)
        assert n_channel_names == n_channels, (
            f'number of channel_names ({n_channel_names}) must match '
            f'number of channels ({n_channels})'
        )
    
    if software == None:
        software = ''
        
    metadata = {
        'Creator': software,
        'Pixels': {
            'PhysicalSizeX': size_x,
            'PhysicalSizeXUnit': unit_x,
            'PhysicalSizeY': size_y,
            'PhysicalSizeYUnit': unit_y,
        },
        'Channel': {'Name': channel_names},
        
    }

    max_size = np.max(img_shape)
    subifds = np.ceil(np.log2(max_size / 1024)).astype(int)
    
    # use optimal tile size for disk space
    tile_size = 16*np.ceil(
        np.array(img_shape[1:]) / (2**subifds) / 16
    ).astype(int)
    options = {
        'tile': tuple(tile_size)
    }

    with tifffile.TiffWriter(output_path, bigtiff=True) as tiff_out:
        tiff_out.write(
            data=out_img,
            metadata=metadata,
            software=software,
            subifds=subifds,
            **options
        )
        for i in range(subifds):
            if i == 0:
                down_2x_img = downsize_img_channels(out_img, is_mask=is_mask)
            else:
                down_2x_img = downsize_img_channels(down_2x_img, is_mask=is_mask)
            tiff_out.write(
                data=down_2x_img,
                subfiletype=1,
                **options
            )

    out_img = out_img.reshape(img_shape_ori)
    return

def downsize_channel(img, is_mask):
    if is_mask:
        return img[::2, ::2]
    else:
        return skimage.transform.downscale_local_mean(img, (2, 2)).astype(img.dtype)

def downsize_img_channels(img, is_mask):
    return np.array([
        downsize_channel(c, is_mask=is_mask)
        for c in img
    ])
