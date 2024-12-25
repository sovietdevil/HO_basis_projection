from ase import Atoms
import numpy as np

def __init__():
    super().__init__()

def generate_column(symbol, thickness, distance, width, gap=0):
    H = distance*(thickness-1)
    #results from pure Aluminum column
    return Atoms(symbol+f'{thickness}', 
               positions=[(width/2, width/2, i*distance + gap) for i in range(thickness)],
               cell=[width, width, H + gap])

def generate_periodic_fcc(symbol, a, nx, ny, nz):
    cell = Atoms(symbol+"4", cell=[a, a, a], pbc=(1, 1, 1), 
            scaled_positions=[(0, 0, 0), (0, 0.5, 0.5), (0.5, 0, 0.5), (0.5, 0.5, 0)])
    supercell = cell * (nx, ny, nz)
    cellprime = cell * (nx+1, ny+1, nz+1)
    cellprime.cell = supercell.cell
    in_cell = [all(position <= [cellprime.cell[i,i] for i in range(3)]) for position in cellprime.positions]
    cell_wrap = cellprime[in_cell]
    return cell_wrap

def set_dopants_column(structure, dopant_list, original_index, changed_index):
    structure_doped = structure.copy()
    structure_doped.set_atomic_number([changed_index if j in dopant_list else original_index for j in range(len(structure))])
    return structure_doped

def set_dopants_fcc(structure, a, dopant_positions, original_index, changed_index):
    dopant_positions = np.array([[a/2, a/2, a*(dopant_positions-1)]])
    dopant_list = np.where(np.isclose(structure.positions, dopant_positions[:,None]).all(-1))[1]
    structure_doped = structure.copy()
    structure_doped.set_atomic_numbers([changed_index if j in dopant_list else original_index for j in range(len(structure))])
    return structure_doped