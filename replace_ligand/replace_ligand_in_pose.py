########################################################################################################################

__doc__ = \
    """
    Main class: ``LigandPoseReplacer``, which adds only an extra method to ``LigandReplacer``.
    But acts as the main conversion method `replace_to_pose`.
    """
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2019 A.D."
__license__ = "MIT"
__version__ = "1"
__citation__ = "N/A"

########################################################################################################################

from .replace_ligand import LigandReplacer
from typing import Optional

import pyrosetta, random

pyrosetta.init('-mute all')


class LigandPoseReplacer(LigandReplacer):
    """
    Can create params and pose it in pyrosetta
    """

    def replace_to_pose(self, probe_name: str, probe_smiles: str, probe_resn: str = 'MMX') -> pyrosetta.Pose:
        pdb = self.replace(probe_smiles=probe_smiles, probe_resn=probe_resn, probe_name=probe_name)
        # pyrosetta.init(f'-extra_res_fa {probe_name}.params -mute all')
        pose = pyrosetta.Pose()
        params_paths = pyrosetta.rosetta.utility.vector1_string()
        params_paths.extend([f"{probe_name}.params"])
        nonstandard_residue_set = pyrosetta.generate_nonstandard_residue_set(pose, params_paths)
        pyrosetta.rosetta.core.import_pose.pose_from_pdbstring(pose, pdb)
        return pose

    def dock(self, pose: pyrosetta.Pose, cst: bool = False) -> None:
        """
        Assumes the chain B is last. Uses ``docking.setup_foldtree`` mysterious method.

        :param pose: replaced structure with ligand
        :return: Changes in place
        """
        dock_jump = 1
        pyrosetta.rosetta.protocols.docking.setup_foldtree(pose, 'A_X', pyrosetta.Vector1([dock_jump]))
        scorefxn = pyrosetta.create_score_function('ligand')
        if cst:
            stm = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
            st = stm.score_type_from_name("atom_pair_constraint")
            scorefxn.set_weight(st, 20)
        docking = pyrosetta.rosetta.protocols.docking.DockMCMProtocol()
        docking.set_scorefxn(scorefxn)
        docking.apply(pose)

    def residue2idx(self, pose, resn):
        res = pyrosetta.rosetta.core.select.residue_selector.ResidueNameSelector(resn)
        # pyrosetta.rosetta.core.select.residue_selector.ChainSelector
        return list(res.apply(pose)).index(1) + 1

    def wiggle(self, pose, idx):
        xyz = pyrosetta.rosetta.numeric.xyzVector_double_t()
        s = 1
        xyz.x = random.gauss(0, s)
        xyz.y = random.gauss(0, s)
        xyz.z = random.gauss(0, s)
        for a in range(1, pose.residue(idx).natoms() + 1):
            pose.residue(idx).set_xyz(a, pose.residue(idx).xyz(a) + xyz)

    def add_cst(self, pose: pyrosetta.Pose,
                res_i: int, atm_i: str, res_j: int, atm_j: str,
                distance: float = 2.4,
                sigma: float = 0.1) \
            -> pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint:
        AtomID = pyrosetta.rosetta.core.id.AtomID
        id_i = AtomID(pose.residue(res_i).atom_index(atm_i), res_i)
        id_j = AtomID(pose.residue(res_j).atom_index(atm_j), res_j)
        ijfunc = pyrosetta.rosetta.core.scoring.constraints.BoundFunc(0.0, distance, sigma, 'cst1')
        cst_ij = pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint(
            id_i, id_j, ijfunc
        )
        pose.add_constraint(cst_ij)
        return cst_ij

    @staticmethod
    def printpair(name: str, apo: pyrosetta.Pose, substrate: pyrosetta.Pose, intermediate: pyrosetta.Pose):
        scorefxn = pyrosetta.get_fa_scorefxn()
        a = scorefxn(apo)
        s = scorefxn(substrate)
        i = scorefxn(intermediate)
        print(f'{name}. Apo {a:.1f}, Substrate {s:.1f} ({s - a:.1f}), intermediate {i:.1f} ({i - a:.1f})')
