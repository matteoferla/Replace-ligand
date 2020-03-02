# replace_ligand package

## Subpackages


* replace_ligand.molfile_to_params package


    * Subpackages


        * replace_ligand.molfile_to_params.rosetta_py package


            * Subpackages


                * replace_ligand.molfile_to_params.rosetta_py.io package


                    * Submodules


                    * replace_ligand.molfile_to_params.rosetta_py.io.mdl_molfile module


                    * replace_ligand.molfile_to_params.rosetta_py.io.pdb module


                    * Module contents


                * replace_ligand.molfile_to_params.rosetta_py.utility package


                    * Submodules


                    * replace_ligand.molfile_to_params.rosetta_py.utility.ForkManager module


                    * replace_ligand.molfile_to_params.rosetta_py.utility.r3 module


                    * replace_ligand.molfile_to_params.rosetta_py.utility.rankorder module


                    * Module contents


            * Module contents


    * Submodules


    * replace_ligand.molfile_to_params.mol_to_params module


    * replace_ligand.molfile_to_params.param_utils module


    * Module contents


## Submodules

## replace_ligand.replace_ligand module

Main class: `LigandReplacer`.
Replace a ligand with another in a PDB from the RCSB PDB. But does solely an RDKit RMSD change of the best conformer.
It is combined with parameterisation.
Gets inherited by `LigandPoseReplacer` to make a pose.
The instatiation is based upon a template, while the replace method gets a probe
and returns a PDB with the ligand best matched.


#### class replace_ligand.replace_ligand.LigandReplacer(code, ref_smiles, ref_resn)
Bases: `object`

Get the PDB block of a conversion:

```python
>>> lr = LigandReplacer(code='3SRG', ref_smiles='O=c1ccc2ccccc2[nH]1', ref_resn='OCH')
>>> print(lr.replace(probe_smiles='Oc1(O)ccc2ccccc2[nH]1', probe_resn='MMX', probe_name='inter'))
```

Requests the resn in the 4


### PDBBlocks()
alias of `pdbblocks`


### \__init__(code, ref_smiles, ref_resn)
The instatiation is based upon a template, while the replace method gets a probe
and returns a PDB with the ligand best matched.


* **Parameters**

    
    * **code** (`str`) – 4 letter PDB code or PDBBlock


    * **ref_smiles** (`str`) – requied form bond order correction (`AllChem.AssignBondOrdersFromTemplate`)


    * **ref_resn** (`str`) – required to extract the ligand



### align_probe_to_target(probe)

* **Parameters**

    **probe** (`Mol`) – modified inplace



* **Return type**

    `int`



* **Returns**

    index of best conformer



### combine(probe_name, probe_resn='MMX')
The probe_resn needs to be untaken by PDB codes for rosetta. see `replace` for more.
Adding -load_PDB_components false during the first initialisation of pyrosetta circumvents this,
but means that native ligands need to declared.


* **Parameters**

    
    * **probe_name** (`str`) – 


    * **probe_resn** (`str`) – 



* **Return type**

    `pdbblocks`



* **Returns**

    


### fill_PDBBlocks(pymol, resn)

* **Return type**

    `pdbblocks`



### get_target_from_pdb()
Fills the `self.target` attribute.


* **Return type**

    `Mol`



* **Returns**

    a `Chem.Mol`



### make_probe(smiles)
This coverts the Smiles string to a RDKit molecule,
for which conformers are generated.
Oddly, using the argument `params=ETKDG()` for
`AllChem.EmbedMultipleConfs` gave results that were problematic,
in that it was very strict in terms of free energy and did not 
consider minor torsions for O-phenylacetate. Whereas openbabel was fine with it.


* **Return type**

    `Mol`



### static mol_with_atom_index(mol)
Returns a copy of the molecule that when displayed shows atom indices.


* **Parameters**

    **mol** (`Mol`) – target molecule



* **Return type**

    `Mol`



* **Returns**

    labelled molecule



### prepare_template()
Part of instatiation for `self.template`.


* **Return type**

    `pdbblocks`



* **Returns**

    a PDBBlocks namedtuple (all, no_ligand, only_ligand).



### replace(probe_smiles, probe_resn, probe_name)
The target probe needs a un-taken PDB 3 letter residue name if the flag `load_PDB_components` is left.
Else the PDB component library version will win!


* **Parameters**

    
    * **probe_smiles** (`str`) – Valid Smiles, protons optional.


    * **probe_resn** (`str`) – untaken 3 letter residue name


    * **probe_name** (`str`) – file friend!



* **Return type**

    `str`



* **Returns**

    

## replace_ligand.replace_ligand_in_pose module

Main class: `LigandPoseReplacer`, which adds only an extra method to `LigandReplacer`.
But acts as the main conversion method replace_to_pose.


#### class replace_ligand.replace_ligand_in_pose.LigandPoseReplacer(code, ref_smiles, ref_resn)
Bases: `replace_ligand.replace_ligand.LigandReplacer`

Can create params and pose it in pyrosetta


### add_cst(pose, res_i, atm_i, res_j, atm_j)

* **Return type**

    `AtomPairConstraint`



### dock(pose, cst)
Assumes the chain B is last. Uses `docking.setup_foldtree` mysterious method.


* **Parameters**

    **pose** (`Pose`) – replaced structure with ligand



* **Return type**

    `None`



* **Returns**

    Changes in place



### static printpair(name, apo, substrate, intermediate)

### replace_to_pose(probe_name, probe_smiles, probe_resn='MMX')

* **Return type**

    `Pose`



### residue2idx(pose, resn)

### wiggle(pose, idx)
## Module contents

See


* `replace_ligand.py` for `LigandReplacer`


* `replace_ligand_in_pose.py` for `LigandPoseReplacer`


* `replace_acetylated.py` for `AcetylPoseReplacer`
