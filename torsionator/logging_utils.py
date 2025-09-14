import logging, sys

def setup_logging(log_file: str):
    logging.basicConfig(
        level=logging.INFO,
        format="%(message)s",
        handlers=[logging.FileHandler(log_file, mode="w"),
                  logging.StreamHandler(sys.stdout)],
        force=True,
    )
    return logging.getLogger("torsionfit")
