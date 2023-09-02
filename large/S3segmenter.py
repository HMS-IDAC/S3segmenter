import argparse 
import sys 
import pathlib 
 
import watershed 
import s3seg_util
import s3seg_qc
import ignored_args

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
        '--probMapChan', 
        default=None,
        type=int,
        nargs="+"
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
    
    ignored_args.add_unsupported_argument(parser)
    args, extra_argv = parser.parse_known_args(argv[1:]) 

    img_path = pathlib.Path(args.imagePath) 
    img_stem = img_path.name.split('.')[0] 
 
    out_path = pathlib.Path(args.outputPath) / img_stem / 'nucleiRing.ome.tif' 
    out_path.parent.mkdir(exist_ok=True, parents=True) 
    
    img_channels = args.probMapChan
    if img_channels is None:
        logging.warning(
            'Image channel used for generating probability maps not specified'
            ' use `--probMapChan CHANNEL` to specify it. Assuming first channel (1)'
        )
        img_channels = [0]
    else: img_channels = [c-1 for c in img_channels]
    assert min(img_channels) >= 0, f'--probMapChannel ({args.probMapChan}) must >= 1'

    if args.pixelSize is not None:
        pixel_size = args.pixelSize
        logging.info(f"Pixel size: {pixel_size} (user supplied)")
    else:
        pixel_size = s3seg_util.detect_pixel_size(img_path)
        if pixel_size is None:
            logging.error(
                'Auto-detect pixel size failed, use `--pixelSize SIZE` to specify it'
            )
            return 1
        logging.info(f"Pixel size: {pixel_size} (from ome-xml)")

    watershed.main([ 
        '', 
        '-i', args.stackProbPath, 
        '-o', str(out_path), 
        '--mcmicro', 
        '--pixel-size', str(pixel_size),
        *extra_argv 
    ]) 

    qc_dir = pathlib.Path(args.outputPath) / img_stem / 'qc' 
    qc_dir.mkdir(exist_ok=True, parents=True)

    s3seg_qc.run_mcmicro(
        out_path,
        qc_dir / f"{img_stem}-nucleiRingOutlines.ome.tif",
        pmap_path=args.stackProbPath,
        img_path=args.imagePath,
        img_channels=img_channels,
        pixel_size=pixel_size
    )
    
    return 0
 
if __name__ == '__main__': 
    sys.exit(main())