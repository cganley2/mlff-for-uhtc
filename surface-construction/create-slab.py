from ase import io
from ase.build import surface, bulk, make_supercell

# locate CIF structure to manipulate
zrb2 = io.read('PATH/TO/REPO/vc-relax-primitive-cell/ZrB2-qe-vc-relax.out')

num_layers = 3 # number of surface layers
vacuum_size = 50 # length of vacuum above surface
n = 4
# Create (001) surface
surface_from_cif = surface(zrb2, (0, 0, 1), layers=num_layers, vacuum=50.0)

# Create supercell
slab = make_supercell(surface_from_cif, [[n, 0, 0], [0, n, 0], [0, 0, 1]])
print('layers: {0}, n={1}; {2} atoms'.format(num_layers, n, len(slab)))

# Uncomment these as needed

# Save as xyz file
# io.write('zrb2-001-slab.xyz'.format(len(slab)), slab, format='xyz')

# Save as LAMMPS data file
# io.write('zrb2-slab-{0}-atoms.data'.format(len(slab)), slab, format='lammps-data')
