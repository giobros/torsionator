import os
from dataclasses import dataclass


@dataclass
class CalculatorFactory:
    obiwan_path: str
    net_charge: int = 0
    spin: int = 1

    def get(self, method: str):
        """Return (ase_calculator, label) for 'mace', 'obi', or 'uma'."""
        if method == "mace":
            from mace.calculators import mace_off
            calc = mace_off(model="medium", device="cuda")
        elif method == "obi":
            from ase.calculators.obi import Obiwan
            model_path = f"{self.obiwan_path}/obiwan_ani1Uani2_FH_VL_2.404"
            calc = Obiwan(model_path=model_path)
        elif method == "uma":
            from fairchem.core import FAIRChemCalculator
            from fairchem.core.units.mlip_unit import load_predict_unit
            checkpoint = os.environ.get("UMA_MODEL_PATH", "uma-s-1p1.pt")
            predictor = load_predict_unit(checkpoint, device="cuda")
            calc = FAIRChemCalculator(predictor, task_name="omol")
            calc._is_uma = True
        else:
            raise ValueError(f"Unknown method: {method!r}")
        calc._net_charge = self.net_charge
        calc._spin = self.spin
        return calc, method


def attach_calc(atoms, calc) -> None:
    """Assign *calc* to *atoms*; inject charge/spin metadata always."""
    atoms.calc = calc
    atoms.info.update({
        "spin": getattr(calc, "_spin", 1),
        "charge": getattr(calc, "_net_charge", 0),
    })
