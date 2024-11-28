import numpy as np
import matplotlib.pyplot as plt
import os
from time import time


class AIAT:
    """
    A class for ab initio atomistic thermodynamics analysis
    and plotting surface energies of materials.
    """

    def __init__(self, energy_dict):
        """
        Initialize the AIAT class with energy parameters.

        Args:
            energy_dict (dict): Contains energies and related properties:
                - 'compound_bulk': Bulk compound energy (eV/formula unit).
                - 'elements_bulk': Bulk energy of individual elements (dict).
                - 'compound_surf': Surface compound energy (list of dicts).
                - 'composition_ratio': Atomic ratios in the compound (dict).
        """
        self.compound_bulk = energy_dict['compound_bulk']
        self.elements_bulk = energy_dict['elements_bulk']
        self.compound_surf = energy_dict['compound_surf']
        self.composition_ratio = energy_dict['composition_ratio']
        self.E_form = self.calculate_formation_energy()

    def calculate_formation_energy(self):
        """
        Calculate the formation energy of the compound.

        Returns:
            float: Formation energy (eV).
        """
        formation_energy = self.compound_bulk
        for element, ratio in self.composition_ratio.items():
            formation_energy -= self.elements_bulk[element] * ratio
        return formation_energy

    def surface_energy_calculator(self, dmu, area, energy, composition):
        """
        Calculate the surface energy for a given surface system.

        Args:
            dmu (array): Chemical potential range (eV).
            area (float): Surface area (Å²).
            energy (float): Surface energy (eV).
            composition (dict): Elemental composition of the surface.

        Returns:
            array: Surface energy values (eV/Å²).
        """
        chemical_potential = sum(
            composition[el] * (dmu if el == 'variable' else self.elements_bulk[el])
            for el in composition
        )
        G = (energy - self.compound_bulk - chemical_potential) / area / 2.0
        return G

    def plot_mu_surf_energy(self, num_points=10, add_mu_range=0.3, ymax=0.3, output_file='aiat.png'):
        """
        Plot the surface energy as a function of chemical potential.

        Args:
            num_points (int): Resolution of the plot.
            add_mu_range (float): Range for chemical potential (eV).
            ymax (float): Maximum y-axis value for the plot.
            output_file (str): Name of the output image file.
        """
        mu_min, mu_max = self.E_form - add_mu_range, add_mu_range
        dmu_range = np.linspace(mu_min, mu_max, num_points)

        print(f"{'Title':<20}\t{'Rich (eV/Å²)':<12}\t{'Lean (eV/Å²)':<12}")
        for system in self.compound_surf:
            title = system['title']
            area = system['area']
            energy = system['energy']
            composition = system['composition']

            G = self.surface_energy_calculator(dmu_range, area, energy, composition)

            plt.plot(dmu_range, G, label=title)
            rich = self.surface_energy_calculator(0, area, energy, composition)[0]
            lean = self.surface_energy_calculator(self.E_form, area, energy, composition)[0]

            print(f"{title:<20}\t{rich:<12.8f}\t{lean:<12.8f}")

        # Plot configurations
        plt.legend()
        plt.plot([0, 0], [-10, 10], '--', color='black')
        plt.plot([self.E_form, self.E_form], [-10, 10], '--', color='black')
        plt.xlim(mu_min, mu_max)
        plt.ylim(0, ymax)
        plt.ylabel('Surface energy (eV/Å²)')
        plt.xlabel('Chemical potential, Δμ (eV)')
        plt.text(0, ymax * 1.01, 'Rich')
        plt.text(self.E_form, ymax * 1.01, 'Lean')
        plt.savefig(output_file)
        print(f"Plot saved to {output_file}")


if __name__ == '__main__':
    # Example input for TiN
    energy_dict = {
        'compound_bulk': -19.628924,  # TiN bulk energy (eV/formula unit)
        'elements_bulk': {
            'Ti': -15.668424 / 2,  # Bulk Ti energy (eV/atom)
            'N': -16.650181 / 2   # Bulk N2 energy per N atom (eV/atom)
        },
        'compound_surf': [
            {'title': 'TiN(100)', 'area': 9.028582, 'energy': -175.35971, 'composition': {'Ti': 9, 'N': 9}},
            {'title': 'TiN(110)', 'area': 25.536674, 'energy': -345.09568, 'composition': {'Ti': 18, 'N': 18}},
            {'title': 'TiN(111):Ti', 'area': 7.8189819, 'energy': -84.675409, 'composition': {'Ti': 5, 'N': 4}},
            {'title': 'TiN(111):N', 'area': 7.8189819, 'energy': -85.44504, 'composition': {'Ti': 4, 'N': 5}},
            {'title': 'TiN(211)', 'area': 22.115394, 'energy': -228.05203, 'composition': {'Ti': 12, 'N': 12}},
            {'title': 'TiN(211)-N', 'area': 22.115394, 'energy': -207.23992, 'composition': {'Ti': 12, 'N': 10}}
        ],
        'composition_ratio': {'Ti': 1, 'N': 1}  # TiN composition
    }

    # Initialize AIAT class
    aiat = AIAT(energy_dict)

    # Timing and execution
    start_time = time()
    aiat.plot_mu_surf_energy()
    end_time = time()

    print(f"\nCode execution completed! Time used: {end_time - start_time:.2f} seconds\n")

    # Open the resulting plot
    os.system(f'xdg-open aiat.png &')
