import argparse 
import sys 
import pathlib 
 
import watershed 
 
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
 
    args, extra_argv = parser.parse_known_args(argv[1:]) 
 
    img_path = pathlib.Path(args.imagePath) 
    img_stem = img_path.name.split('.')[0] 
 
    out_path = pathlib.Path(args.outputPath) / img_stem / 'nuclei.ome.tif' 
    out_path.parent.mkdir(exist_ok=True, parents=True) 
 
    return watershed.main([ 
        '', 
        '-i', args.stackProbPath, 
        '-o', str(out_path), 
        '--mcmicro', 
        *extra_argv 
    ]) 

 
if __name__ == '__main__': 
    sys.exit(main())