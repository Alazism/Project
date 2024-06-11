# -*- coding: utf-8 -*-
"""
Created on Fri May 31 16:41:12 2024

@author: kumul
"""

import math

def calculate_properties(Mass_fraction_core, Delta_ppm_shift, Delta_ppm_error, Electron_effective_mass, core_density,
                         chi_vol_bulk, MW_ligand, shell_density, chi_mol_ligand, Mass_fraction_ligand,
                         MW_solvent, Density_solvent, chi_mol_solvent, Bohrs, free_space):
    
    # Constants
    pi = math.pi
    
    # Step 1: chi_mass_obs calculation
    chi_mass_obs = (3 * Delta_ppm_shift) / (4 * pi * 10**6)
    chi_mass_obs_error = (3 * Delta_ppm_error) / (4 * pi * 10**6)
    
    # Step 2: chi_vol_obs calculation
    Mass_fraction_shell = 1 - Mass_fraction_core
    chi_vol_obs = ((Mass_fraction_core / core_density) + (Mass_fraction_shell / shell_density))**(-1) * chi_mass_obs
    chi_vol_obs_error = chi_vol_obs * (chi_mass_obs_error / chi_mass_obs)
    
    # Step 3: mass_lig_core calculation
    mass_lig_core = Mass_fraction_ligand / Mass_fraction_core
    
    # Step 4: chi_vol_core calculation
    chi_vol_core = (1 + (core_density / shell_density) * mass_lig_core)**(-1) * chi_vol_bulk
    
    # Step 5: chi_vol_shell calculation
    chi_vol_shell = (1 + (shell_density / core_density) * mass_lig_core**-1)**(-1) * (chi_mol_ligand * shell_density) / MW_ligand
    
    # Step 6: chi_vol_solvent calculation
    chi_vol_solvent = (chi_mol_solvent * Density_solvent) / MW_solvent
    
    # Step 7: chi_vol_conduction calculation
    chi_vol_conduction = chi_vol_obs - chi_vol_core - chi_vol_shell + chi_vol_solvent
    chi_vol_conduction_error = chi_vol_obs_error  # Errors for chi_vol_core, chi_vol_shell, and chi_vol_solvent are assumed to be negligible
    
    # Step 8: chi_vol_Pauli calculation
    chi_vol_Pauli = chi_vol_conduction / (1 - 1/3 * Electron_effective_mass**-2)
    chi_vol_Pauli_error = chi_vol_conduction_error / (1 - 1/3 * Electron_effective_mass**-2)
    
    # Step 9: chi_vol_landau calculation
    chi_vol_landau = -(1 / Electron_effective_mass)**2 * chi_vol_Pauli / 3
    chi_vol_landau_error = abs(chi_vol_landau * (chi_vol_Pauli_error / chi_vol_Pauli))
    
    # Step 10: gE_f calculation
    gE_f = chi_vol_Pauli / (free_space * Bohrs**2)
    gE_f_error = abs(gE_f * (chi_vol_Pauli_error / chi_vol_Pauli))
    
    # Print all calculated values with errors
    print(f"chi_mass_obs: {chi_mass_obs:.6e} ± {chi_mass_obs_error:.6e}")
    print(f"chi_vol_obs: {chi_vol_obs:.6e} ± {chi_vol_obs_error:.6e}")
    print(f"mass_lig_core: {mass_lig_core:.6e}")
    print(f"chi_vol_core: {chi_vol_core:.6e}")
    print(f"chi_vol_shell: {chi_vol_shell:.6e}")
    print(f"chi_vol_solvent: {chi_vol_solvent:.6e}")
    print(f"chi_vol_conduction: {chi_vol_conduction:.6e} ± {chi_vol_conduction_error:.6e}")
    print(f"chi_vol_Pauli: {chi_vol_Pauli:.6e} ± {chi_vol_Pauli_error:.6e}")
    print(f"chi_vol_landau: {chi_vol_landau:.6e} ± {chi_vol_landau_error:.6e}")
    print(f"gE_f: {gE_f:.6e} ± {gE_f_error:.6e}")

# input values (CHANGE THESE FOR SYSTEM)
Mass_fraction_core = 0.477950848
Delta_ppm_shift = 1.512896831
Delta_ppm_error = 0.018858825

# Change for metal
Electron_effective_mass = 2.2
core_density = 12.023
chi_vol_bulk = -3.00e-5

# Change for ligand
MW_ligand = 201.392
shell_density = 0.845
chi_mol_ligand = -0.00016025


Mass_fraction_ligand = 1- Mass_fraction_core

# Change for solvent
MW_solvent = 120.38
Density_solvent = 1.5
chi_mol_solvent = -0.00005993

# Constants
Bohrs = 9.27401e-24
free_space = 0.000001256636

calculate_properties(Mass_fraction_core, Delta_ppm_shift, Delta_ppm_error, Electron_effective_mass, core_density,
                     chi_vol_bulk, MW_ligand, shell_density, chi_mol_ligand, Mass_fraction_ligand,
                     MW_solvent, Density_solvent, chi_mol_solvent, Bohrs, free_space)
