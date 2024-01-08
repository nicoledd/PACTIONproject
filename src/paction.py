import sys
import argparse

from algo import runAlgo
from GTtoData import getGTdata
from GTtoInput import getGTinput
from compare import runCompare


def main(args):
    GTdata = getGTdata(args.g)
    GTinput = getGTinput(args.g)
    Tdata = runAlgo(GTinput)
    runCompare(GTdata, Tdata, args.o)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-o', type=str, help='output prefix', required=True)
    parser.add_argument('-g', type=str, help='ground truth', required=True)

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    main(args)

