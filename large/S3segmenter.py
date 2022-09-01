import argparse 
import sys 
import pathlib 
 
import watershed 
import s3seg_util
import logging
 
def main(argv=sys.argv): 
 
    parser = argparse.ArgumentParser( 
        description=( 
            's3segmentor for large images' 
        ), 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter 
    ) 
 
    parser.add_argument( 
        '--imagePath', 
        default=argparse.SUPPRESS, 
        required=True 
    ) 
    parser.add_argument( 
        '--stackProbPath', 
        default=argparse.SUPPRESS, 
        required=True 
    ) 
    parser.add_argument( 
        '--outputPath', 
        default='.' 
    ) 
    parser.add_argument( 
        '--pixelSize', 
        default=None,
        type=float 
    )
 
    args, extra_argv = parser.parse_known_args(argv[1:]) 
 
    img_path = pathlib.Path(args.imagePath) 
    img_stem = img_path.name.split('.')[0] 
 
    out_path = pathlib.Path(args.outputPath) / img_stem / 'nuclei.ome.tif' 
    out_path.parent.mkdir(exist_ok=True, parents=True) 
 
    return watershed.main([ 
    if args.pixelSize is not None:
        pixel_size = args.pixelSize
    else:
        pixel_size = s3seg_util.detect_pixel_size(img_path)
        if pixel_size is None:
            logging.error(
                'Auto-detect pixel size failed, use `--pixelSize SIZE` to specify it'
            )
            return 1

    watershed.main([ 
        '', 
        '-i', args.stackProbPath, 
        '-o', str(out_path), 
        '--mcmicro', 
        '--pixel-size', str(pixel_size),
        *extra_argv 
    ]) 

 
if __name__ == '__main__': 
    sys.exit(main())