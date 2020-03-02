########################################################################################################################

__doc__ = \
    """
    Main class: ``LigandReplacer``.
    Replace a ligand with another in a PDB from the RCSB PDB. But does solely an RDKit RMSD change of the best conformer.
    It is combined with parameterisation.
    Gets inherited by ``LigandPoseReplacer`` to make a pose.
    The instatiation is based upon a template, while the `replace` method gets a probe
    and returns a PDB with the ligand best matched.
    
    """
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2019 A.D."
__license__ = "MIT"
__version__ = "1"
__citation__ = "N/A"

########################################################################################################################

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdFMCS, rdmolfiles
import pymol2
from . import molfile_to_params
from collections import namedtuple
from typing import Tuple, Sequence


class LigandReplacer:
    """
    Get the PDB block of a conversion:

    >>> lr = LigandReplacer(code='3SRG', ref_smiles='O=c1ccc2ccccc2[nH]1', ref_resn='OCH')
    >>> print(lr.replace(probe_smiles='Oc1(O)ccc2ccccc2[nH]1', probe_resn='MMX', probe_name='inter'))

    Requests the resn in the 4
    """

    PDBBlocks = namedtuple('pdbblocks', ['all', 'no_ligand', 'only_ligand'])

    def __init__(self, code: str, ref_smiles: str, ref_resn: str):
        """
        The instatiation is based upon a template, while the `replace` method gets a probe
        and returns a PDB with the ligand best matched.

        :param code: 4 letter PDB code or PDBBlock
        :param ref_smiles: requied form bond order correction (``AllChem.AssignBondOrdersFromTemplate``)
        :param ref_resn: required to extract the ligand
        """
        self.code = code
        self.ref_resn = ref_resn
        self.ref_smiles = ref_smiles
        # create self.template: namedtuple('template', ['all', 'no_ligand', 'only_ligand']) of pdb blocks
        self.template = self.prepare_template()
        # creates self.target: Chem.Mol corrected based on self.ref_smiles
        self.target = self.get_target_from_pdb()

    def prepare_template(self) -> PDBBlocks:
        """
        Part of instatiation for ``self.template``.

        :return: a PDBBlocks namedtuple (all, no_ligand, only_ligand).
        """
        with pymol2.PyMOL() as pymol:
            if len(self.code) == 4:
                pymol.cmd.fetch(self.code)
            elif '.pdb' in self.code:  # file.
                pymol.cmd.load(self.code)
            else:
                raise ValueError
            pymol.cmd.remove('solvent')
            pymol.cmd.h_add('*')
            pymol.cmd.alter('all','segi=""')
            pymol.cmd.sort()
            return self.fill_PDBBlocks(pymol, self.ref_resn)

    def fill_PDBBlocks(self, pymol: pymol2.PyMOL, resn: str) -> PDBBlocks:
        # pymol is assumed to be started and the data loaded.
        return self.PDBBlocks(all=pymol.cmd.get_pdbstr('*'),
                              no_ligand=pymol.cmd.get_pdbstr(f'not resn {resn}'),
                              only_ligand=pymol.cmd.get_pdbstr(f'resn {resn}'))

    def get_target_from_pdb(self) -> Chem.Mol:
        """
        Fills the ``self.target`` attribute.

        :return: a ``Chem.Mol``
        """
        # The reference molecule
        ref = Chem.MolFromSmiles(self.ref_smiles)
        # The PDB conformations
        target = Chem.MolFromPDBBlock(self.template.only_ligand)  # , sanitize=False, strictParsing=False)
        target = AllChem.AssignBondOrdersFromTemplate(ref, target)
        # target = Chem.AddHs(target) #done at the PDB step
        return target

    def replace(self, probe_smiles: str, probe_resn: str, probe_name: str) -> str:
        """
        The target probe needs a un-taken PDB 3 letter residue name if the flag ``load_PDB_components`` is left.
        Else the PDB component library version will win!

        :param probe_smiles: Valid Smiles, protons optional.
        :param probe_resn: untaken 3 letter residue name
        :param probe_name: file friend!
        :return:
        """
        probe = self.make_probe(probe_smiles)
        cid = self.align_probe_to_target(probe)
        aligned_file = f'{probe_name}.aligned.mol'
        #Chem.MolToMolFile(probe, aligned_file, confId=cid)
        writer = rdmolfiles.SDWriter(aligned_file)
        writer.write(probe, confId=cid)
        for i in range(probe.GetNumConformers()):
            if i == cid:
                continue
            writer.write(probe, confId=i)
        writer.close()
        molfile_to_params.run(aligned_file,
                              conformers_in_one_file=True,
                              name=probe_resn,
                              amino_acid=None,
                              chain='X',
                              pdb=probe_name,
                              clobber=True)
        pdbblocks = self.combine(probe_name, probe_resn)
        return pdbblocks.all

    def make_probe(self, smiles: str) -> Chem.Mol:
        """
        This coverts the Smiles string to a RDKit molecule,
        for which conformers are generated.
        Oddly, using the argument ``params=ETKDG()`` for
        ``AllChem.EmbedMultipleConfs`` gave results that were problematic,
        in that it was very strict in terms of free energy and did not 
        consider minor torsions for O-phenylacetate. Whereas openbabel was fine with it.

        .. code-block:: bash
             obabel -:"CC(=O)Oc1ccccc1" --gen3d -p -osdf -O temp.sdf
             obabel temp.sdf -osdf -O conf.sdf --confab --conf 100 --original --rcutoff 0.05

        """
        probe = Chem.MolFromSmiles(smiles)
        assert probe is not None, f'{smiles} is dodgy.'
        probe = Chem.AddHs(probe)
        #         AllChem.EmbedMolecule(probe)
        #         AllChem.UFFOptimizeMolecule(probe, maxIters=2000)

        #         #Chem.rdPartialCharges.ComputeGasteigerCharges(probe)
        idx = AllChem.EmbedMultipleConfs(probe,
                                         numConfs=100,
                                         forceTol=0.2,
                                         pruneRmsThresh=0.1)
        #         for i in idx:
        #             AllChem.UFFOptimizeMolecule(probe, maxIters=2000, confId=i)
        AllChem.UFFOptimizeMoleculeConfs(probe)
        rdMolAlign.AlignMolConformers(probe)
        assert len(idx) > 0, 'No conformers generated!'
        return probe

    def align_probe_to_target(self, probe: Chem.Mol) -> int:  # implace
        """

        :param probe: modified inplace
        :return: index of best conformer
        """
        ### find what is common
        common = self._get_common(probe)
        ### Align them
        overlap_target = self.target.GetSubstructMatch(common)
        overlap_probe = probe.GetSubstructMatch(common)
        atomMap = [(probe_at, target_at) for probe_at, target_at in zip(overlap_probe, overlap_target)]
        rmss = [rdMolAlign.AlignMol(probe,
                                    self.target,
                                    prbCid=i,
                                    atomMap=atomMap,
                                    maxIters=500) for i in range(probe.GetNumConformers())]
        # print(rmss)
        best_i = rmss.index(min(rmss))
        return best_i

    def _get_common(self, probe: Chem.Mol) -> Chem.Mol:
        res = Chem.rdFMCS.FindMCS([probe, self.target]  # ,
                                  # matchValences=True, threshold=0.1
                                  # atomCompare=Chem.rdFMCS.AtomCompare.CompareElements #,
                                  # bondCompare=Chem.rdFMCS.BondCompare.CompareAny #CompareOrder
                                  )
        return Chem.MolFromSmarts(res.smartsString)

    def combine(self, probe_name: str, probe_resn: str = 'MMX') -> PDBBlocks:
        """
        The probe_resn needs to be untaken by PDB codes for rosetta. see ``replace`` for more.
        Adding -load_PDB_components false during the first initialisation of pyrosetta circumvents this,
        but means that native ligands need to declared.

        :param probe_name:
        :param probe_resn:
        :return:
        """
        with pymol2.PyMOL() as pymol:
            pymol.cmd.read_pdbstr(self.template.no_ligand, 'nolig')
            pymol.cmd.load(f'{probe_name}.pdb')
            return self.fill_PDBBlocks(pymol, probe_resn)

    @staticmethod
    def mol_with_atom_index(mol: Chem.Mol) -> Chem.Mol:
        """
        Returns a copy of the molecule that when displayed shows atom indices.

        :param mol: target molecule
        :return: labelled molecule
        """
        cp = Chem.Mol(mol)
        atoms = cp.GetNumAtoms()
        for idx in range(atoms):
            cp.GetAtomWithIdx(idx).SetProp('molAtomMapNumber', str(mol.GetAtomWithIdx(idx).GetIdx()))
        return cp
