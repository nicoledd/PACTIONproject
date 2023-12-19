import sys
import argparse

from progressive_caller import solveProgressivePaction
from simultaneous_caller import solveSimultaneousPaction
from segment_caller import solveSegmentPaction


# python3 paction.py -m progressive -p test_suite/example_input/props_mode0.out test_suite/example_input/props_mode1.out test_suite/example_input/props_mode2.out -t test_suite/example_input/tree_mode0.tsv None None -o test_suite/example_output/


def main(args):
    assert (args.m=='progressive' or args.m == 'simultaneous' or args.m == 'both' or args.m == 'segment'), 'not a valid algorithm'
    if args.m == 'progressive':
        solveProgressivePaction(args)
    elif args.m == 'simultaneous':
        solveSimultaneousPaction(args)
    elif args.m == 'both':
        solveProgressivePaction(args)
        solveSimultaneousPaction(args)
    elif args.m == 'segment':
        solveSegmentPaction(args)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    # assume first modality is snv (props+tree), the rest are copy-number states (no trees)
    parser.add_argument('-m', type=str, help='which algorithm to run? [\'progressive\',\'simultaneous\',\'both\']', required=True)
    parser.add_argument('-c', type=str, help='cna files', required=True, nargs='*')
    parser.add_argument('-s', type=str, help='snv file', required=True)
    parser.add_argument('-o', type=str, help='output prefix', required=True)
    parser.add_argument('-g', type=str, help='ground truth', required=False)


    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    main(args)

