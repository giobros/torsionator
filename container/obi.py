import os
import numpy as np
import tensorflow as tf
from ase.calculators.calculator import Calculator
from ase.io import read
import sys

obi_root = os.getenv("OBI_ROOT")
if obi_root and obi_root not in sys.path:
    sys.path.insert(0, obi_root)

from architectures import net_utils
eh_to_eV = 27.2114
class Obiwan(Calculator):
    implemented_properties = ['energy', 'forces']  

    def __init__(self, model_path, model_name="obiwan", **kwargs):
        """
        Inizializza il calcolatore Obiwan.

        Parameters:
        - model_path: str, percorso ai pesi del modello.
        - model_name: str, nome del modello (default: 'obiwan').
        """
        super().__init__(**kwargs)
        self.model = net_utils.getModel(model_name=model_name)
        self.model.loadWeights(model_path, dynamic_mode=False)

    def calculate(self, atoms=None, properties=['energy'], system_changes=None):
        """
        Calcola energia e forze utilizzando il modello Obiwan.
        """
        super().calculate(atoms, properties, system_changes)

      
        species = atoms.get_chemical_symbols()
        coordinates = atoms.get_positions()

    
        species_tensor = tf.constant(species, dtype=tf.string)
        coordinates_tensor = tf.constant(coordinates, dtype=tf.float32)

      
        batched_coordinates = tf.expand_dims(coordinates_tensor, axis=0)
        batched_species = tf.expand_dims(species_tensor, axis=0)

     
        energy, forces = self.model.computeEnergyAndForces([batched_coordinates, batched_species])


        self.results['energy'] = eh_to_eV * energy.numpy().item()
        self.results['forces'] = eh_to_eV * forces.numpy().squeeze()

       
