from dataclasses import dataclass
from typing import Tuple
from mace.calculators import mace_off
from ase.calculators.obi import Obiwan

@dataclass
class CalculatorFactory:
    obiwan_path: str

    def get(self, method: str):
        """Return (ase_calculator, label) for 'mace' | 'obi'."""
        if method == "mace":
            return mace_off(model="medium", device="cpu"), "mace"
        elif method == "obi":
            model_path = f"{self.obiwan_path}/obiwan_ani1Uani2_FH_VL_2.404"
            return Obiwan(model_path=model_path), "obi"
        raise ValueError(f"Unknown method: {method}")

    def label(self, method: str) -> str:
        """Return the label for the given method."""
        mapping = {"mace": "mace", "obi": "obi"}
        try:
            return mapping[method]
        except KeyError:
            raise ValueError(f"Unknown method: {method}")
