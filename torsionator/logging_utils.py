import logging
import sys
import warnings


_SILENT_LOGGERS = [
    "ase",
    "mace",
    "torch",
    "torch.distributed",
    "faiss",
    "urllib3",
    "requests",
    "numba",
    "e3nn",
]


def setup_logging(log_file: str) -> logging.Logger:
    logging.basicConfig(
        level=logging.INFO,
        format="%(message)s",
        handlers=[
            logging.FileHandler(log_file, mode="w"),
            logging.StreamHandler(sys.stdout),
        ],
        force=True,
    )

    for name in _SILENT_LOGGERS:
        logging.getLogger(name).setLevel(logging.ERROR)

    warnings.filterwarnings("ignore")

    return logging.getLogger("torsionfit")
