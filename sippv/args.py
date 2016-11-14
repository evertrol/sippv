import argparse


def parse():
    def nargs_range(nmin, nmax):
        class NargsRange(argparse.Action):
            def __call__(self, parser, args, values, option_string=None):
                if not nmin <= len(values) <= nmax:
                    msg = '"{f}" requires between {nmin} and {nmax} arguments'.format(
                        f=self.dest, nmin=nmin, nmax=nmax)
                    raise argparse.ArgumentTypeError(msg)
                setattr(args, self.dest, values)
        return NargsRange

    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+',
                        help="Input FITS files")
    parser.add_argument('-o', '--order', type=int,
                        help="PV order")
    parser.add_argument('-r', '--radial', action='store_true',
                        help="Use the radial term")
    parser.add_argument('-n', '--ndata', type=int, default=[64, 64],
                        action=nargs_range(1, 2), nargs='+',
                        help="Use n x n pixels to solve")
    parser.add_argument('--ext', action='append',
                        help="Extension to correct; "
                        "option can be used multiple times. "
                        "Default is the primary hdu.")
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help="Verbose level")
    args = parser.parse_args()

    if len(args.ndata) == 1:
        args.ndata = [args.ndata[0], args.ndata[0]]
    if not args.ext:
        args.ext = [0]

    return args
