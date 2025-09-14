# unit conversions
CONV_EV_TO_EH = 0.0367493
CONV_EH_TO_KCAL_MOL = 627.509474

# Lennard-Jones parameters and vdw radii
LJ_PARAMETERS = {
    'H': {'sigma': 2.5,  'epsilon': 0.044},
    'C': {'sigma': 3.4,  'epsilon': 0.359},
    'N': {'sigma': 3.25, 'epsilon': 0.71},
    'O': {'sigma': 2.96, 'epsilon': 0.88},
    'S': {'sigma': 3.56, 'epsilon': 1.046},
    'P': {'sigma': 3.74, 'epsilon': 0.669},
    'F': {'sigma': 2.94, 'epsilon': 0.61},
    'Cl': {'sigma': 3.52, 'epsilon': 1.108},
    'Br': {'sigma': 3.73, 'epsilon': 1.519},
    'I': {'sigma': 4.009, 'epsilon': 2.098},
}
VDW_RADII = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52,
    "F": 1.47, "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98,
}