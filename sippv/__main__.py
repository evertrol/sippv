from . import core
from .args import parse as parse_args
from .logging import setup as setup_logger


if __name__ == '__main__':
    import sys
    import os
    if sys.argv[0] == __file__:
        sys.argv[0] = "{executable} -m sippv".format(executable=sys.executable)

    args = parse_args()
    setup_logger('sippv', args.verbose)
    core.run(args.files, order=args.order, radial=args.radial,
             ndata=args.ndata, extensions=args.ext)
