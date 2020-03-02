# Replace-ligand
Example of how to replace a ligand with a similar in Rdkit and Pyrosetta.

## Atomic name replacer
In a seperate repo, [github.com/matteoferla/AtomicRenamer](https://github.com/matteoferla/AtomicRenamer),
I have put a script that given a SMILES string, creates a mol2 file,
where the atom names match closely to that of a template.
This useful for downstream steps, like using pyrosetta where the atom names are important,
and one wants them consistent without bending over backwards.

## Ligand replacer
> This requies my modded [mol_to_params.py](https://github.com/matteoferla/mol_to_params.py).
> mol_to_params is to make a Rosetta Params file for the novel ligand.
> This version is modded, namely, a Python3 port with a change that allows it to be called easily like a normal module.
> Unfortunately due to licence I cannot share the  publically,

The class `LigandReplacer` give a PDB filename or code, 
a reference smiles (to correct the bond order) of the ligand in the template complex
and the 3 letter name of the ligand, returns an object where the method fix, give a probe name (for files),
an untaken 3 letter name and the smiles of the ligand will replace it to the other by generating and 
choosing the best conformer (which is saved as first conformer in the SDF also).

    lr = LigandReplacer(code='template.pdb', ref_smiles='O=c1ccc2ccccc2[nH]1', ref_resn='OCH')
    dihydroxycoumarine = lr.replace(probe_name='dihydroxycoumarine',
                                    probe_resn='WQW',
                                    probe_smiles='O=c1ccc2ccccc2o1')


The 3 letter code is a bit of a pain.
It has to be untaken if PDB components are loaded (needed for ions, ATP etc.).
Try 3 random letters: https://www.rcsb.org/ligand/WQU.
An improvement could be adding in the saving of the multiconformer SDF a step to save only conformers with an 
RMSD from previously saved ones greater than 0.1. Also the conformer generation is slightly odd,
namely using:

    AllChem.EmbedMultipleConfs(probe,
                               numConfs=100,
                               forceTol=0.2,
                               pruneRmsThresh=0.1)
    AllChem.UFFOptimizeMoleculeConfs(probe)

As opposed to using `params=AllChem.ETKDG()`.
The reason for this is that I wanted to add some torsion in the ligand, even if it makes it unhappy.

The second class, `LigandPoseReplacer`, inherits `LigandReplacer` and adds a few extra.
Instead of `.replace(..)`, calling `.replace_to_pose(probe_name, probe_smiles, probe_resn)`
will return a pyrosetta pose.
This pose can be docked.
Beforehand `.add_cst(pose, pose_idx_A, atom_name_A, pose_idx_B, atom_name_B, distance, sigma)`
can add a harmonic atom pair constraint to the pose.
`pose_idx` is the pose number.
To convert from PDB to pose index in pyrosetta you use `pdb2pose = pose.pdb_info().pdb2pose` as `pdb2pose(res:int=n, chain:str=c)`.
Then `.dock(pose, cst=True)` can be run.

I won't lie. For most applications, the replacement method 
`align_probe_to_target(self, probe: Chem.Mol) -> int` will need to be overridden by a copy, most often 
changing the atom mapping (atomMap) of `rdMolAlign.AlignMol` and adding a `weights` list after doing some weird replacement.