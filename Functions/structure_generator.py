from ase import Atoms
from ase.build import stack
from matplotlib import pyplot as plt
from ase.visualize import view
import numpy as np

def construction_periodic(structure, nx, ny, nz):
    '''
    Construct periodic structure with boundary included
    input:
    structure: ASE Atoms structure within a unit cell;
    nx, ny, nz: numbers of repetitions in x, y, z directions
    return:
    cell_wrap: the ASE Atoms structure of the goal supercell
    '''
    supercell = structure * (nx, ny, nz)
    cellprime = structure * (nx+1, ny+1, nz+1)
    cellprime.cell = supercell.cell
    in_cell = [all(position <= [cellprime.cell[i,i] for i in range(3)]) for position in cellprime.positions]
    cell_wrap = cellprime[in_cell]
    return cell_wrap

def set_dopants(structure, dopant_positions, original_index, changed_index):
    '''
    Generate structure according to a certain index or a list of indexes
    '''
    dopant_list = np.where(np.isclose(structure.positions, dopant_positions[:,None]).all(-1))[1]
    structure_doped = structure.copy()
    structure_doped.set_atomic_numbers([changed_index if j in dopant_list else original_index for j in range(len(structure))])
    return structure_doped

def construct_amorphous(element, number, a, b, c):
    coordx = np.random.random(number)
    coordy = np.random.random(number)
    coordz = np.random.random(number)
    coords = np.transpose(np.vstack([coordx, coordy, coordz]))
    return Atoms(element+f"{number}", 
                 scaled_positions=coords,
                 cell=[a, b, c])