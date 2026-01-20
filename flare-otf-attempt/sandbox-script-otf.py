# This code is based on the tutorial at: https://github.com/mir-group/flare/blob/master/tutorials/sparse_gp_tutorial.ipynb

# Import numpy and matplotlib
import numpy as np
from numpy.random import random
import matplotlib.pyplot as plt
import matplotlib

# Increase the matplotlib font size.
font = {'size': 22}

matplotlib.rc('font', **font)

# flare++ imports
from flare.bffs.sgp import SGP_Wrapper
from flare.bffs.sgp.calculator import SGP_Calculator
from flare.bffs.sgp._C_flare import B2, NormalizedDotProduct, SparseGP, Structure

# flare imports
from flare.learners.otf import OTF
from flare.io import otf_parser

# ASE imports
import ase
from ase import Atoms, units
from ase.calculators.espresso import Espresso
from ase.build import surface, make_supercell
from ase.visualize import view
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, \
    Stationary, ZeroRotation
from ase.build import fcc111, add_adsorbate
from ase.io import write
from ase.io import read


################### EXAMPLE 2 - otf

# set calculator with settings
input_data = {
    'pseudo_dir': '/home/cganley2/qe-data/pseudopotentials/SSSP_1.3.0_PBE_precision',
    'tprnfor': True,
    'tstress': True,
    'calculation': 'scf',
    'occupations': 'smearing',
    'smearing': 'gaussian',
    'degauss': 0.01,
    'ecutwfc': 50,
    'ecutrho': 400,
    'conv_thr': 1e-7
}
kpts = (2, 2, 1)

pseudopotentials = {'Zr': 'Zr.pbe-spn-kjpaw_psl.1.0.0.UPF', 'B': 'B.pbe-n-kjpaw_psl.1.0.0.UPF'}
print('line 47')
calc = Espresso(command='mpirun -np 4 pw.x -npool 1 -in espresso.pwi > espresso.pwo', input_data=input_data, pseudopotentials=pseudopotentials, kpts=kpts)

# choose the initial structure, arbitrarily the midpoint of AIMD
# Load training data
# zrb2 = ase.io.read('/home/cganley2/qe-data/zrb2/structures/ZrB2_conventional_standard_mp-1472.cif')
print('line 52')
# # Create (001) surface
# surface_from_cif = surface(zrb2, (0, 0, 1), layers=3, vacuum=10.0, periodic=True)
# # Create supercell
# n = 3
# atoms = make_supercell(surface_from_cif, [[n, 0, 0], [0, n, 0], [0, 0, 1]])
atoms = ase.io.read('/home/cganley2/scitech2026/pre-training-data/zrb2-slab-001/qe-calcs/relax-9atoms/zrb2-slab-001-9atoms-relax.out')
n_atoms = len(atoms)

# choose MD settings
# Set MD parameters.
md_engine = "VelocityVerlet"
md_dict = {}

# Set the initial velocity to 300 K.
temperature = 500  # in K
MaxwellBoltzmannDistribution(atoms, temperature_K=temperature)
Stationary(atoms)  # zero linear momentum
ZeroRotation(atoms)  # zero angular momentum

print('line 71')
# Choose model settings
# Create sparse GP model.
species_map = {40: 0, 5: 1}  # Aluminum (atomic number 13) is species 0
cutoff = 5.0  # in A
sigma = 2.0  # in eV
power = 2  # power of the dot product kernel
kernel = NormalizedDotProduct(sigma, power)
cutoff_function = "quadratic"
many_body_cutoffs = [cutoff]
radial_basis = "chebyshev"
radial_hyps = [0., cutoff]
cutoff_hyps = []
n_species = 2
N = 8
lmax = 3
descriptor_settings = [n_species, N, lmax]
descriptor_calculator = B2(
  radial_basis,
  cutoff_function,
  radial_hyps,
  cutoff_hyps,
  descriptor_settings
)

print('line 96')
# Set the noise values.
sigma_e = 0.001 * n_atoms  # eV (1 meV/atom)
sigma_f = 0.05  # eV/A
sigma_s = 0.0006  # eV/A^3 (about 0.1 GPa)

# Choose uncertainty type.
# Other options are "DTC" (Deterministic Training Conditional) or
# "SOR" (Subset of Regressors).
variance_type = "local"  # Compute uncertainties on local energies (normalized)

# Choose settings for hyperparameter optimization.
max_iterations = 20  # Max number of BFGS iterations during optimization
opt_method = "L-BFGS-B"  # Method used for hyperparameter optimization

# Bounds for hyperparameter optimization.
# Keeps the energy noise from going to zero.
bounds = [(None, None), (sigma_e, None), (None, None), (None, None)]

# Create a model wrapper that is compatible with the flare code.
gp_model = SGP_Wrapper(
    [kernel],
    [descriptor_calculator],
    cutoff,
    sigma_e,
    sigma_f,
    sigma_s,
    species_map,
    variance_type=variance_type,
    energy_training=True,
    force_training=True,
    stress_training=False,
    opt_method=opt_method,
    bounds=bounds,
    max_iterations=max_iterations,
)

# Create an ASE calculator based on the GP model.
flare_calculator = SGP_Calculator(gp_model)

print('line 134')
# Choose OTF settings
# Set up OTF object.
init_atoms = list(range(n_atoms))  # Initial environments to include in the sparse set
output_name = 'zrb2'  # Name of the output file
std_tolerance_factor = -0.01  # Uncertainty tolerance for calling QM
train_hyps = (0, 10)  # Freeze hyperparameter optimization after second QM call
min_steps_with_model = 0  # Min number of steps between DFT calls
update_style = "threshold"  # Strategy for adding sparse environments
update_threshold = 0.005  # Threshold for determining which sparse environments to add
force_only = False  # Train only on forces or include energies and stresses
rescale_steps = [1000, 2000, 3000, 4000]
rescale_temps = [1500, 2500, 3000, 4000]

otf_params = {
    'init_atoms': init_atoms,
    'output_name': output_name,
    'std_tolerance_factor': std_tolerance_factor,
    'train_hyps': train_hyps,
    'min_steps_with_model': min_steps_with_model,
    'update_style': update_style,
    'update_threshold': update_threshold,
    'rescale_steps': rescale_steps,
    'rescale_temps': rescale_temps,
}

# Create OTF object.
timestep = 0.001  # units of ps
number_of_steps = 5000
test_otf = OTF(
    atoms,
    timestep,
    number_of_steps,
    calc,
    md_engine,
    md_dict,
    flare_calc=flare_calculator,
    force_only=force_only,
    **otf_params,
)

print('line 175')
# Run on-the-fly dynamics.
test_otf.run()

print('line 179')
# Parse the output file.
output_file = '{0}.out'.format(output_name)
otf_trajectory = otf_parser.OtfAnalysis(output_file)

# Plot temperature and energy vs. simulation time.
times = otf_trajectory.times
eam_times = otf_trajectory.dft_times

temps = otf_trajectory.thermostat['temperature']
eam_temps = otf_trajectory.gp_thermostat['temperature']

gp_energies = otf_trajectory.thermostat['potential energy']
eam_energies = otf_trajectory.gp_thermostat['potential energy']

fig, ax = plt.subplots((1, 2))

ax[0].plot(times, temps)
ax[0].plot(eam_times, eam_temps, 'kx')
ax[0].xlabel('Time (ps)')
ax[0].ylabel('Temperature (K)')

ax[1].plot(times, gp_energies)
ax[1].plot(eam_times, eam_energies, 'kx')
ax[1].xlabel("Time (ps)")
ax[1].ylabel("Potential energy (eV)")

fig.savefig('{0}-plots.png'.format(output_name), dpi=400)
#################### EXAMPLE 1 - aspirin

# # Load training data
# sic_bulk_aimd = ase.io.read('/home/cganley2/qe-data/sic-3c/simulations/aimd-2500K/SiC_3C-aimd-2500K.out', index=':')
# n_strucs = len(sic_bulk_aimd)
# forces = np.array([value.get_forces() for index, value in enumerate(sic_bulk_aimd)])
# positions = np.array([value.get_positions() for index, value in enumerate(sic_bulk_aimd)])
# noa = len(sic_bulk_aimd[0].get_chemical_symbols())
# cell = np.array(sic_bulk_aimd[0].get_cell())

# species = sic_bulk_aimd[0].get_chemical_symbols()
# species_code = {'Si': 0, 'C': 1}

# coded_species = []
# for spec in species:
#     coded_species.append(species_code[str(spec)])

# # Choose training and validation structures.
# training_size = 100
# validation_size = 20
# np.random.seed(1)
# shuffled_frames = [int(n) for n in range(n_strucs)]
# np.random.shuffle(shuffled_frames)

# training_pts = shuffled_frames[0:training_size]
# validation_pts = shuffled_frames[training_size:training_size+validation_size]

# # Specify the hyperparameters of the Sparse GP model

# # Define many-body descriptor.
# cutoff = 3.7  # A
# n_species = 2  # Silicon, Carbon
# N = 12  # Number of radial basis functions
# lmax = 3  # Largest L included in spherical harmonics
# radial_basis = "chebyshev"  # Radial basis set
# cutoff_name = "quadratic"  # Cutoff function
# radial_hyps = [0, cutoff]
# cutoff_hyps = []
# descriptor_settings = [n_species, N, lmax]

# # Define a B2 object.
# B2_descriptor = B2(radial_basis, cutoff_name, radial_hyps, cutoff_hyps,
#                  descriptor_settings)

# # The GP class can take a list of descriptors as input, but here
# # we'll use a single descriptor.
# descriptors = [B2_descriptor]

# # Define kernel function as normalized dot product kernel
# sigma = 2.0
# power = 2
# dot_product_kernel = NormalizedDotProduct(sigma, power)

# # Define a list of kernels.
# # There needs to be one kernel for each descriptor.
# kernels = [dot_product_kernel]

# # Define sparse GP. These noise values are initialized to the expected error level for each quantity
# sigma_e = 0.12 * noa  # Energy noise (in kcal/mol, so about 5 meV/atom)
# sigma_f = 0.115  # Force noise (in kcal/mol/A, so about 5 meV/A)
# sigma_s = 0.014  # Stress noise (in kcal/A^3, so about 0.1 GPa)
# gp_model = SparseGP(kernels, sigma_e, sigma_f, sigma_s)

# # Calculate descriptors of the validation and training structures.
# print("Computing descriptors of validation points...")
# validation_strucs = []
# validation_forces = np.zeros((validation_size, noa, 3))
# for n, snapshot in enumerate(validation_pts):
#     pos = positions[snapshot]
#     frcs = forces[snapshot]

#     # Create structure object, which computes and stores descriptors.
#     struc = \
#         Structure(cell, coded_species, pos, cutoff, descriptors)
#     validation_strucs.append(struc)
#     validation_forces[n] = frcs
# print("Done.")

# print("Computing descriptors of training points...")
# training_strucs = []
# training_forces = np.zeros((training_size, noa, 3))
# for n, snapshot in enumerate(training_pts):
#     pos = positions[snapshot]
#     frcs = forces[snapshot]

#     # Create structure object, which computes and stores descriptors.
#     struc = \
#         Structure(cell, coded_species, pos, cutoff, descriptors)

#     # Assign force labels to the training structure.
#     struc.forces = frcs.reshape(-1)

#     training_strucs.append(struc)
#     training_forces[n] = frcs
# print("Done.")



# # Train the model.
# print("Training the GP...")
# batch_size = 10  # monitor the MAE after adding this many frames
# n_strucs = np.zeros(batch_size)
# mb_maes = np.zeros(batch_size)
# for m in range(training_size):
#     train_struc = training_strucs[m]

#     # Add training structure and sparse environments.
#     gp_model.add_training_structure(train_struc)
#     gp_model.add_all_environments(train_struc)

#     if (m + 1) % batch_size == 0:
#         # Update the sparse GP training coefficients.
#         gp_model.update_matrices_QR()

#         # Predict on the validation set.
#         pred_forces = np.zeros((validation_size, noa, 3))
#         for n, test_struc in enumerate(validation_strucs):
#             gp_model.predict_SOR(test_struc)
#             pred_vals = test_struc.mean_efs[1:-6].reshape(noa, 3)
#             pred_forces[n] = pred_vals

#         # Calculate and store the MAE.
#         batch_no = int((m + 1) / batch_size)
#         mae = np.mean(np.abs(validation_forces - pred_forces))
#         n_strucs[batch_no - 1] = batch_size * batch_no
#         mb_maes[batch_no - 1] = mae
#         print("Batch %i MAE: %.2f kcal/mol/A" % (batch_no, mae))

##################### END EXAMPLE 1 - aspirin