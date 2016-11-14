import logging


def setup(name=None, verbose=0):
    logger = logging.getLogger(name)
    formatter = logging.Formatter("%(message)s")
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    handler.setLevel('DEBUG')
    logger.addHandler(handler)
    logger.setLevel(['WARNING', 'INFO', 'DEBUG'][verbose])
