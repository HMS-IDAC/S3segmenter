import tifffile 
import numpy as np 
import cv2 
import zarr 
import palom
 
import skimage.morphology 
import skimage.util 
import skimage.measure 
import skimage.segmentation 
 
import dask_image.ndmeasure 
 
import dask.array as da 
import dask.diagnostics 
 
import pathlib 
import logging 
import pprint 
logging.basicConfig( 
    format="%(asctime)s | %(levelname)-8s | %(message)s (%(filename)s:%(lineno)s)", 
    datefmt="%Y-%m-%d %H:%M:%S", 
    level=logging.INFO 
)

import _version
 
 
def filter_label_area(label_img, area_min, area_max): 
    if np.all(label_img == 0): 
        return label_img 
    ids, counts = np.unique(label_img, return_counts=True) 
    filtered = np.where( 
        (counts >= area_min) & (counts < area_max), 
        ids, 
        0 
    ) 
    indexer = np.arange(ids.max()+1) 
    indexer[ids] = filtered 
    return indexer[label_img].astype(np.int32) 
 
 
def filter_label_intensity(label_img, intensity_img, intensity_min): 
    if np.all(label_img == 0): 
        return label_img 
    table = skimage.measure.regionprops_table( 
        label_img, 
        intensity_img, 
        properties=('label', 'mean_intensity') 
    ) 
    ids, ints = table['label'], table['mean_intensity'] 
    filtered = np.where( 
        ints >= intensity_min, 
        ids, 
        0 
    ) 
    indexer = np.arange(ids.max()+1) 
    indexer[ids] = filtered 
    return indexer[label_img].astype(np.int32) 
 
 
def gaussian_2d_cv2(img, sigma): 
    return cv2.GaussianBlur( 
        img, (0, 0), sigma,  
        borderType=cv2.BORDER_REFLECT 
    ) if sigma != 0 else img 
 
 
 
class WatershedSegmentor: 
 
    def __init__( 
        self, 
        contour_img, 
        intensity_img, 
        h=0.01, 
        sigma=1, 
        footprint_size=3, 
        overlap_size=None, 
        mean_intensity_min=0, 
        area_min=0, 
        area_max=np.inf,
        pixel_size=None
    ) -> None: 
        self.contour_img = contour_img 
        self._contour_img = contour_img 
        self.intensity_img = intensity_img 
        self._intensity_img = intensity_img 
        self.h = h 
        self.sigma = sigma 
        self.footprint_size = footprint_size 
        self.overlap_size = overlap_size 
        self.mean_intensity_min = mean_intensity_min 
        self.area_min = area_min 
        self.area_max = area_max 
        self.pixel_size = pixel_size
 
        self.mask = None 
        self.config = None 
        self.tests = [] 
        self._update_config() 
 
    def run(self, config_id=None, compute=True): 
        self._update_config() 
        if config_id is not None: self.use_config(config_id) 
        logging.info('Run initiating') 
        logging.info( 
            'Configuration:\n\n{}\n\n'.format( 
                pprint.pformat(self.config[max(self.config.keys())]) 
            ) 
        ) 
        filtered = self.filter_intensity() 
        logging.info('Start configuring `label`') 
        re_labeled = dask_image.ndmeasure.label(filtered)[0] 
        logging.info('End configuring `label`') 
        if compute: 
            with dask.diagnostics.ProgressBar(): 
                return re_labeled.persist() 
        return re_labeled 
 
    def write(self, file_path, img=None, config_id=None): 
        file_path = pathlib.Path(file_path) 
        file_name = file_path.name 
        if img is None: img = self.run(config_id, compute=False) 
        if file_name.endswith('.zarr'): 
            logging.info(f'Writing to {file_path}') 
            with dask.diagnostics.ProgressBar(): 
                return img.to_zarr(file_path) 
        if file_name.endswith(('.ome.tiff', '.ome.tif')): 
            logging.info(f'Writing to {file_path}')
            pixel_size = self.pixel_size
            if self.pixel_size is None:
                pixel_size = 1
            return palom.pyramid.write_pyramid(
                [img],
                file_path,
                pixel_size=pixel_size,
                downscale_factor=2,
                compression='zlib',
                is_mask=True,
                tile_size=1024,
                save_RAM=True,
                kwargs_tifffile=dict(software=f"s3segmenter-large v{_version.VERSION}")
            )
        logging.warning('Write failed: output file type not supported') 
        return 
     
    def write_expanded(self, file_path, expand_size, img=None): 
        file_path = pathlib.Path(file_path) 
        expand_size = int(expand_size) 
        if isinstance(img, (str, pathlib.Path)): 
            in_path = pathlib.Path(img) 
            in_name = in_path.name 
            if in_name.endswith('.zarr'): 
                import zarr 
                img = da.from_zarr(zarr.open(in_path)) 
            elif in_name.endswith(('.ome.tiff', '.ome.tif')): 
                # maybe should use `chunks=2048`? 
                img = da.from_array(tifffile.imread(in_path), chunks=2048) 
            else: 
                logging.warning(f'Image reading ({in_path}) failed. No file written to {file_path}') 
                return 
        expanded = img.map_overlap( 
            skimage.segmentation.expand_labels, 
            distance=expand_size, 
            depth=self.overlap_size, 
            boundary='none', 
            dtype=np.int32 
        ) 
        return self.write(file_path, img=expanded) 
 
    def gaussian_img(self, img=None, sigma=None): 
        if img is None: 
            img = self.contour_img 
            if np.issubdtype(img.dtype, np.integer): 
                img = np.invert(img) 
            else: 
                img *= -1 
        if sigma is None: 
            sigma = self.sigma 
        fimg = img.map_blocks(skimage.util.img_as_float32, dtype=np.float32) 
        if sigma == 0: 
            return fimg 
        return fimg.map_overlap( 
            gaussian_2d_cv2, 
            sigma=sigma, 
            depth=4*int(sigma), 
            dtype=np.float32, 
            boundary='none' 
        ) 
 
    def peak_img(self, img=None, h=None, footprint_size=None, overlap_size=None): 
        if img is None: img = self.gaussian_img() 
        if h is None: h = self.h 
        if footprint_size is None: footprint_size = self.footprint_size 
        if overlap_size is None: overlap_size = self.overlap_size 
        return img.map_overlap( 
            skimage.morphology.extrema.h_maxima, 
            h=h, 
            footprint=np.ones((footprint_size, footprint_size)), 
            depth=overlap_size, 
            dtype=bool, 
            boundary='none' 
        ) > 0 
         
    def segmented_img(self, img=None, seed_img=None, mask=None, overlap_size=None): 
        if img is None: img = self.contour_img 
        if seed_img is None:  
            seed_img_binary = self.peak_img() 
            logging.info('Start configuring `label`') 
            seed_img = dask_image.ndmeasure.label(seed_img_binary)[0] 
            logging.info('End configuring `label`') 
        if mask is None: mask = self.mask 
        if overlap_size is None: overlap_size = self.overlap_size 
        return da.map_overlap( 
            skimage.segmentation.watershed, 
            img, 
            seed_img, 
            watershed_line=True, 
            mask=mask, 
            depth=overlap_size, 
            boundary='none', 
            dtype=np.int32 
        ) 
     
    def filter_area(self, label_img=None, area_min=None, area_max=None, overlap_size=None): 
        if label_img is None: label_img = self.segmented_img() 
        if area_min is None: area_min = self.area_min 
        if area_max is None: area_max = self.area_max 
        if overlap_size is None: overlap_size = self.overlap_size 
        if (area_min == 0) and (area_max == np.inf): 
            return label_img 
        return da.map_overlap( 
            filter_label_area, 
            label_img, 
            area_min=area_min, 
            area_max=area_max, 
            depth=overlap_size, 
            boundary='none', 
            dtype=np.int32 
        ) 
 
    def filter_intensity( 
        self, label_img=None, intensity_img=None, 
        intensity_min=None, overlap_size=None 
    ): 
        if label_img is None: label_img = self.filter_area() 
        if intensity_img is None: intensity_img = self.intensity_img 
        if intensity_min is None: intensity_min = self.mean_intensity_min 
        if overlap_size is None: overlap_size = self.overlap_size 
        return da.map_overlap( 
            filter_label_intensity, 
            label_img, 
            intensity_img, 
            intensity_min=intensity_min, 
            depth=overlap_size, 
            boundary='none', 
            dtype=np.int32 
        ) 
     
    def _update_config(self): 
        config = dict( 
            h=self.h, 
            sigma=self.sigma, 
            footprint_size=self.footprint_size, 
            overlap_size=self.overlap_size, 
            mean_intensity_min=self.mean_intensity_min, 
            area_min=self.area_min, 
            area_max=self.area_max 
        ) 
        if self.config is None: 
            self.config = {0: config} 
        else: 
            last_key = max(self.config.keys()) 
            last = self.config[last_key] 
            if last != config: 
                self.config.update({last_key+1: {**last, **config}}) 
 
    def test_img(self, funcs, contour_img=None, intensity_img=None): 
        if contour_img is None: contour_img = self._contour_img 
        if intensity_img is None: intensity_img = self._intensity_img 
        self.contour_img = contour_img 
        self.intensity_img = intensity_img 
        try: 
            out = None 
            for f in funcs: 
                out = f() 
        except Exception as e: 
            print('Test failed') 
            print(e) 
            out = None 
        finally: 
            self.contour_img = self._contour_img 
            self.intensity_img = self._intensity_img 
        self._update_config() 
        self.tests.append( 
            (self.config[max(self.config.keys())], out) 
        ) 
        return out 
     
    def use_config(self, config_id): 
        config = self.config[config_id] 
        for k in config: 
            setattr(self, k, config[k]) 
        self._update_config() 
 
 
 
import argparse 
import pathlib 
import sys 
import os 
import dask 
import gc 
import dask.system 
 
 
def main(argv=sys.argv): 
 
    parser = argparse.ArgumentParser( 
        description=( 
            'Watershed segment probability maps' 
        ), 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter 
    )  
    parser.add_argument( 
        '-i', 
        metavar='input-file', 
        help='path of input probability map', 
        default=argparse.SUPPRESS, 
        required=True 
    ) 
    parser.add_argument( 
        '-o', 
        metavar='output-ome-path', 
        help='full path to the output ome-tiff file', 
        default=argparse.SUPPRESS, 
        required=True 
    ) 
    parser.add_argument( 
        '--gaussian-sigma', 
        help='sigma value to gaussian blur probability map for peak finding', 
        default=1.0,
        type=float
    ) 
    parser.add_argument( 
        '--maxima-h', 
        help='`h` in `skimage.morphology.extrema.h_maxima`', 
        default=0.01, 
        type=float 
    ) 
    parser.add_argument( 
        '--maxima-footprint-size', 
        help='size of the square fooptrint used in `skimage.morphology.extrema.h_maxima`', 
        default=3, 
        type=int 
    ) 
    parser.add_argument( 
        '--overlap-size', 
        help='pixels overlapping for block computation', 
        default=128, 
        type=int 
    ) 
    parser.add_argument( 
        '--mean-intensity-min', 
        help='minimal mean intensity in each segmment', 
        default=128, 
        type=float 
    ) 
    parser.add_argument( 
        '--area-min', 
        help='minimal num of pixels (inclusive) in each segmment', 
        default=10, 
        type=int 
    ) 
    parser.add_argument( 
        '--area-max', 
        help='maximal num of pixels (exclusive) in each segmment', 
        default=500, 
        type=int 
    ) 
    parser.add_argument( 
        '--expand-size', 
        help='pixels to expand from each segment', 
        default=0, 
        type=int 
    ) 
    parser.add_argument( 
        '--mcmicro', 
        help='flag to name output files in mcmicro format', 
        default=False, 
        action='store_true' 
    ) 
    parser.add_argument( 
        '--pixel-size', 
        help='image pixel size in micron', 
        default=None, 
        type=float
    ) 
     
    args = parser.parse_args(argv[1:]) 
     
    if len(argv) == 1: 
        parser.print_help() 
        return 0 
 
    assert args.o.endswith(('.zarr', '.ome.tif', '.ome.tiff')) 
 
    print() 
    logging.info(f'Reading {args.i}') 
    probability_maps = palom.reader.OmePyramidReader(args.i).pyramid[0][:2].compute()
    logging.info(f'Probability map shape: {probability_maps.shape}') 

    pixel_size = args.pixel_size
    if pixel_size is None:
        try:
            pixel_size = palom.reader.OmePyramidReader(args.i).pixel_size
        except Exception as err:
            print(err)
            pixel_size = 1.0
            logging.warning(
                f"Pixel size not specified, using {pixel_size} µm as a placeholder"
            ) 
 
    segmentor = WatershedSegmentor( 
        da.from_array(probability_maps[1], chunks=2048), 
        da.from_array(probability_maps[0], chunks=2048), 
        args.maxima_h, 
        args.gaussian_sigma, 
        args.maxima_footprint_size, 
        args.overlap_size, 
        args.mean_intensity_min, 
        args.area_min, 
        args.area_max,
        args.pixel_size
    ) 
 
    segmentor.write(args.o) 
 
    if args.expand_size != 0: 
        probability_maps = None 
        segmentor = None 
        gc.collect() 
 
        logging.info(f'Expanding {args.expand_size} pixels') 
        path_expanded = None 
        if args.mcmicro: 
            path_o = pathlib.Path(args.o) 
            path_expanded = path_o.parent / 'cellRing.ome.tif' 
        path_expanded = expand_mask_from_file( 
            args.o, 
            args.expand_size, 
            output_path=path_expanded 
        ) 
 
        logging.info(f'Generating difference mask') 
        path_difference = None 
        if args.mcmicro: 
            path_o = pathlib.Path(args.o) 
            path_difference = path_o.parent / 'cytoRing.ome.tif' 
        difference_mask_from_file( 
            args.o, 
            path_expanded, 
            output_path=path_difference 
        ) 
    
    logging.info('Done') 
    return 0 
 
def expand_mask_from_file(
    input_path, expand_size, output_path=None, pixel_size=None
): 
    expanded_path = output_path 
    if expanded_path is None: 
        input_path = pathlib.Path(input_path) 
        suffix = ''.join(input_path.suffixes) 
        expanded_path = input_path.parent / input_path.name.replace( 
            suffix, f'-expanded_{expand_size}{suffix}' 
        ) 
    segmentor = WatershedSegmentor( 
        None, None, overlap_size=2*expand_size, pixel_size=pixel_size
    ) 
    z = zarr.open( 
        tifffile.imread(input_path, aszarr=True, series=0, level=0) 
    ) 
     
    segmentor.write_expanded( 
        expanded_path, 
        expand_size, 
        da.from_zarr(z).rechunk(2048) 
    ) 
    return expanded_path 
 
def difference_mask_from_file(
    path_1, path_2, output_path=None, pixel_size=None
): 
    if output_path is None: 
        input_path = pathlib.Path(path_1) 
        suffix = ''.join(input_path.suffixes) 
        stem_2 = pathlib.Path(path_2).name.split('.')[0] 
        output_path = input_path.parent / input_path.name.replace( 
            suffix, f'-{stem_2}{suffix}' 
        ) 
    z1, z2 = [ 
        zarr.open( 
            tifffile.imread(p, aszarr=True, series=0, level=0) 
        ) 
        for p in (path_1, path_2) 
    ] 
 
    out_mask = da.map_blocks( 
        # this assumes one mask is a subset of the other mask; computing 
        # `np.sum(x>0)` instead of `x.sum()` to prevent int overflow 
        lambda x, y: (x != y)*x if np.sum(x>0) > np.sum(y>0) else (x != y)*y, 
        da.from_zarr(z1), 
        da.from_zarr(z2) 
    ) 
    segmentor = WatershedSegmentor( 
        None, None, pixel_size=pixel_size
    ) 
    segmentor.write( 
        file_path=output_path, 
        img=out_mask 
    ) 
    return 
 
 
if __name__ == '__main__': 
    num_workers = dask.system.cpu_count() 
    dask.config.set(scheduler='threads', num_workers=num_workers) 
    sys.exit(main())
