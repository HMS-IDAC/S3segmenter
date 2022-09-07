import argparse


def get_old_parser():
    parser = argparse.ArgumentParser( 
        description=( 
            's3segmentor for large images' 
        ), 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter 
    ) 
    parser.add_argument("--imagePath")
    parser.add_argument("--contoursClassProbPath",default ='')
    parser.add_argument("--nucleiClassProbPath",default ='')
    parser.add_argument("--stackProbPath",default ='')
    parser.add_argument("--outputPath")
    parser.add_argument("--dearrayPath")
    parser.add_argument("--maskPath")
    parser.add_argument("--probMapChan",type = int, default = -1)
    parser.add_argument("--mask",choices=['TMA', 'tissue','none'],default = 'tissue')
    parser.add_argument("--crop",choices=['interactiveCrop','autoCrop','noCrop','dearray','plate'], default = 'noCrop')
    parser.add_argument("--cytoMethod",choices=['hybrid','distanceTransform','bwdistanceTransform','ring'],default = 'distanceTransform')
    parser.add_argument("--nucleiFilter",choices=['IntPM','LoG','Int','none'],default = 'IntPM')
    parser.add_argument("--nucleiRegion",choices=['watershedContourDist','watershedContourInt','watershedBWDist','dilation','localThreshold','localMax','bypass','pixellevel'], default = 'watershedContourInt')
    parser.add_argument("--pixelThreshold",type = float, default = -1)
    parser.add_argument("--segmentCytoplasm",choices = ['segmentCytoplasm','ignoreCytoplasm'],default = 'ignoreCytoplasm')
    parser.add_argument("--cytoDilation",type = int, default = 5)
    parser.add_argument("--logSigma",type = int, nargs = '+', default = [3, 60])
    parser.add_argument("--CytoMaskChan",type=int, nargs = '+', default=[2])
    parser.add_argument("--pixelMaskChan",type=int, nargs = '+', default=[2])
    parser.add_argument("--TissueMaskChan",type=int, nargs = '+', default=0)
    parser.add_argument("--detectPuncta",type=int, nargs = '+', default=[0])
    parser.add_argument("--punctaSigma", nargs = '+', type=float, default=[0])
    parser.add_argument("--punctaSD", nargs = '+', type=float, default=[4])
    parser.add_argument("--saveMask",action='store_false')
    parser.add_argument("--saveFig",action='store_false')
    return parser


def add_unsupported_argument(curr_parser, old_parser=get_old_parser()):

    old_options = set([
        k 
        for k, v in old_parser._option_string_actions.items() 
        if not issubclass(type(v), argparse._HelpAction)
    ])
    curr_options = set([
        k 
        for k, v in curr_parser._option_string_actions.items() 
        if not issubclass(type(v), argparse._HelpAction)
    ])
    ignored_options = old_options - curr_options
    for o in ignored_options:
        curr_parser.add_argument(
            o, action=IgnoreAction, default=argparse.SUPPRESS, nargs='*', metavar=''
        )
    
    mark_ignored_help_strings(curr_parser)
    return


# https://gist.github.com/bsolomon1124/44f77ed2f15062c614ef6e102bc683a5
import logging


class IgnoreAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        for s in self.option_strings:
            logging.warning(f"Argument {s} is not supported and is ignored")


def mark_ignored_help_strings(parser, prefix="IGNORED"):
    for action in parser._actions:
        if isinstance(action, IgnoreAction):
            h = action.help
            if h is None:
                action.help = prefix
            else:
                action.help = prefix + ": " + h

