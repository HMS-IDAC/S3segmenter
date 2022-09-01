import ome_types

def detect_pixel_size(img_path):
    try:
        metadata = ome_types.from_tiff(img_path)
        pixel_size = metadata.images[0].pixels.physical_size_x
    except Exception as err:
        print(err)
        print()
        print('Pixel size detection using ome-types failed')
        pixel_size = None
    return pixel_size