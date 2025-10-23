import os
from dataclasses import dataclass

@dataclass
class Config:
    base_dir: str = "/data"
    obiwan_path: str = os.getenv("OBIWAN_MODEL_PATH", "")
    log_file: str = None

    def __post_init__(self):
        if self.log_file is None:
            self.log_file = os.path.join(self.base_dir, "workflow.log")