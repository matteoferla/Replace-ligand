########################################################################################################################

__doc__ = \
    """
    Modified mol_to_params
    """
__author__ = "Modified by Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2019 A.D."
__license__ = "Rosetta's"
__version__ = "2"
__citation__ = "N/A"

########################################################################################################################


from types import SimpleNamespace
from typing import Optional, Union
from .mol_to_params import core
from optparse import OptionParser, IndentedHelpFormatter

import os, sys, copy, io, random #{{{
if not hasattr(sys, "version_info") or sys.version_info < (3,5):
    raise ValueError("Script requires Python 3.5 or higher!")

def run(infile,
        name='LIG',
        pdb: str = 'ligand',
        centroid=False,
        chain: Optional[str] = 'X',
        center:str=None,
        max_confs: Optional[int] = 5000,
        root_atom: Optional[int] = None,
        nbr_atom=None,
        kinemage=None,
        amino_acid: Union[None, bool] = None,
        clobber: bool = False,
        no_param: bool = False,
        no_pdb: bool = False,
        extra_torsion_output: bool = False,
        keep_names: bool = False,
        long_names: bool = False,
        recharge: Union[None, bool] = None,
        m_ctrl=None,
        mm_as_virt: bool = False,
        skip_bad_conformers: bool = False,
        conformers_in_one_file: bool = False
        ):
    """
    infile: INPUT.mol | INPUT.sdf | INPUT.mol2
    name: name ligand residues NM1,NM2,... instead of LG1,LG2,...
    pdb: prefix for PDB file names
    centroid: write files for Rosetta centroid mode too
    chain: The chain letter to use for the output PDB ligand.
    center: translate output PDB coords to have given heavy-atom centroid ('X,Y,Z')
    max_confs: don't expand proton chis if above this many total confs
    root_atom: which atom in the molfile is the root? (indexed from 1)
    nbr_atom: which atom in the molfile is the nbr atom? (indexed from 1)
    kinemage: write ligand topology to FILE
    amino_acid: set up params file for modified amino acid; .mol2 only; edit chis afterward.  Implies --keep-names.
    clobber: overwrite existing files
    no_param: skip writing .params files (for debugging)
    no_pdb: skip writing .pdb files (for debugging)
    extra_torsion_output: writing additional torsion files
    keep_names: leaves atom names untouched except for duplications
    long_names: if specified name is longer than 3 letters, keep entire name in param NAME field instead of truncating
    recharge: ignore existing partial charges, setting total charge to CHG
    m_ctrl: read additional M control lines from FILE
    mm_as_virt: assign mm atom types as VIRT, rather than X
    skip_bad_conformers: If a conformer has atoms in the wrong order, skip it and continue rather than dying
    conformers_in_one_file: Output 1st conformer to NAME.pdb and all others to NAME_conformers.pdb

    Original documentation.
    Converts a small molecule in an MDL Molfile with "M SPLT" and "M ROOT"
    records into a series of .params residue definition files for Rosetta.
    Also writes out the ligand conformation as PDB HETATMs.

    If an SD file is given as input instead, the first entry is used for
    generating topology / parameter files, and they all are used for
    generating PDB-style coordinates in separate, numbered files.
    These multiple files can optionally be concatenated into a single file,
    which can then be specified with an additional PDB_ROTAMERS line in the
    .params file to include the extra conformations as ligand rotamers.
    Multiple models may also be supplied in MOL2 format, which does not support
    M ROOT and M SPLT records but does allow for partial charges.
    File type is deduced from the extension.

    To divide a ligand into fragments by "breaking" bonds (optional):
    M SPLT atom_no1 atom_no2

    To specify a neighbor atom for a ligand fragment (optional):
    M NBR atom_no

    To specify a root atom for a ligand fragment (optional):
    M ROOT atom_no

    The "M" records (M SPLT, M NBR, M ROOT) can alternatively be specified in
    a separate control file, which can be used with MOL2 format files.

    Note that for ligands with multiple rotamers, Rosetta overlays the ligands
    based on the neighbor atom (not the root atom), such that the position of the
    neighbor atom and the orientation of the atoms bonded to the neighbor atom is
    the same. When using ligand rotamers, it is recommended to confirm that the
    neighbor atom falls in an appropriate position.

    Expects that the input ligand has already had aromaticity "perceived",
    i.e. that it contains aromatic bonds rather than alternating single and double
    bonds (Kekule structure).

    Optionally writes a kinemage graphics visualization of the atom tree,
    neighbor atom selection, fragments, etc -- very helpful for debugging
    and for visualizing exactly what was done to the ligand.
    """
    fields = ['name', 'pdb', 'centroid', 'chain', 'center', 'max_confs', 'root_atom', 'nbr_atom', 'kinemage',
              'amino_acid', 'clobber', 'no_param', 'no_pdb', 'extra_torsion_output', 'keep_names', 'long_names',
              'recharge', 'm_ctrl', 'mm_as_virt', 'skip_bad_conformers', 'conformers_in_one_file']

    options = SimpleNamespace(**{k: v for k, v in locals().items() if k in fields})
    # namedtuple does not work as it has to change.
    core(infile, options)

# Better handle multiple paragraph descriptions.
class PreformattedDescFormatter (IndentedHelpFormatter):
    def format_description(self, description):
        return description.strip() + "\n" # Remove leading/trailing whitespace

def main(argv): #{{{
    """
Converts a small molecule in an MDL Molfile with "M SPLT" and "M ROOT"
records into a series of .params residue definition files for Rosetta.
Also writes out the ligand conformation as PDB HETATMs.

If an SD file is given as input instead, the first entry is used for
generating topology / parameter files, and they all are used for
generating PDB-style coordinates in separate, numbered files.
These multiple files can optionally be concatenated into a single file,
which can then be specified with an additional PDB_ROTAMERS line in the
.params file to include the extra conformations as ligand rotamers.
Multiple models may also be supplied in MOL2 format, which does not support
M ROOT and M SPLT records but does allow for partial charges.
File type is deduced from the extension.

To divide a ligand into fragments by "breaking" bonds (optional):
M SPLT atom_no1 atom_no2

To specify a neighbor atom for a ligand fragment (optional):
M NBR atom_no

To specify a root atom for a ligand fragment (optional):
M ROOT atom_no

The "M" records (M SPLT, M NBR, M ROOT) can alternatively be specified in
a separate control file, which can be used with MOL2 format files.

Note that for ligands with multiple rotamers, Rosetta overlays the ligands
based on the neighbor atom (not the root atom), such that the position of the
neighbor atom and the orientation of the atoms bonded to the neighbor atom is
the same. When using ligand rotamers, it is recommended to confirm that the
neighbor atom falls in an appropriate position.

Expects that the input ligand has already had aromaticity "perceived",
i.e. that it contains aromatic bonds rather than alternating single and double
bonds (Kekule structure).

Optionally writes a kinemage graphics visualization of the atom tree,
neighbor atom selection, fragments, etc -- very helpful for debugging
and for visualizing exactly what was done to the ligand.
    """ # Preformatted
    parser = OptionParser(usage="usage: %prog [flags] { INPUT.mol | INPUT.sdf | INPUT.mol2 }", formatter=PreformattedDescFormatter())
    parser.set_description(main.__doc__)
    # parser.add_option("-short", ["--long"],
    #   action="store|store_true|store_false",
    #   default=True|False|...
    #   type="string|int|float",
    #   dest="opt_name",
    #   help="store value in PLACE",
    #   metavar="PLACE",
    # )
    parser.add_option("-n", "--name",
        default="LG",
        help="name ligand residues NM1,NM2,... instead of LG1,LG2,...",
        metavar="NM"
    )
    parser.add_option("-p", "--pdb",
        default=None, # same as --name, see below
        help="prefix for PDB file names",
        metavar="FILE"
    )
    parser.add_option("-c", "--centroid",
        default=False,
        action="store_true",
        help="write files for Rosetta centroid mode too"
    )
    parser.add_option("--chain",
        default='X',
        type="string",
        help="The chain letter to use for the output PDB ligand."
    )
    parser.add_option("--center",
        default=None, # same as --name, see below
        help="translate output PDB coords to have given heavy-atom centroid",
        metavar="X,Y,Z"
    )
    parser.add_option("-m", "--max-confs",
        default=5000, # 400 (default Omega max) * 9 (one sp3 H with -ex1) = 3600
        type="int",
        help="don't expand proton chis if above this many total confs",
        metavar="MAX"
    )
    parser.add_option("--root_atom",
        type="int",
        help="which atom in the molfile is the root? (indexed from 1)",
        metavar="ATOM_NUM"
    )
    parser.add_option("--nbr_atom",
        type="int",
        help="which atom in the molfile is the nbr atom? (indexed from 1)",
        metavar="ATOM_NUM"
    )
    parser.add_option("-k", "--kinemage",
        default=None,
        help="write ligand topology to FILE",
        metavar="FILE"
    )
    parser.add_option("-a", "--amino-acid",
        default=None,
        help="set up params file for modified amino acid; .mol2 only; edit chis afterward.  Implies --keep-names.",
        metavar="ALA"
    )
    parser.add_option("--clobber",
        default=False,
        action="store_true",
        help="overwrite existing files"
    )
    parser.add_option("--no-param",
        default=False,
        action="store_true",
        help="skip writing .params files (for debugging)"
    )
    parser.add_option("--no-pdb",
        default=False,
        action="store_true",
        help="skip writing .pdb files (for debugging)"
    )
    parser.add_option("--extra_torsion_output",
        default=False,
        action="store_true",
        help="writing additional torsion files"
    )
    parser.add_option("--keep-names",
        default=False,
        action="store_true",
        help="leaves atom names untouched except for duplications"
    )
    parser.add_option("--long-names",
        default=False,
        action="store_true",
        help="if specified name is longer than 3 letters, keep entire name in param NAME field instead of truncating"
    )
    parser.add_option("--recharge",
        type="int",
        help="ignore existing partial charges, setting total charge to CHG",
        metavar="CHG",
    )
    parser.add_option("--m-ctrl",
        default=None,
        help="read additional M control lines from FILE",
        metavar="FILE"
    )
    parser.add_option("--mm-as-virt",
        default=False,
        dest="mm_as_virt",
        help="assign mm atom types as VIRT, rather than X",
        action="store_true"
    )
    parser.add_option("--skip-bad-conformers",
        default=False,
        dest="skip_bad_conformers",
        help="If a conformer has atoms in the wrong order, skip it and continue rather than dying",
        action="store_true"
    )
    parser.add_option("--conformers-in-one-file",
        default=False,
        help="Output 1st conformer to NAME.pdb and all others to NAME_conformers.pdb",
        action="store_true"
    )


    (options, args) = parser.parse_args(args=argv)
    if options.pdb is None: options.pdb = options.name
    if options.amino_acid is not None: options.keep_names = True

    if len(args) < 1:
        parser.print_help()
        print("Must specify input .mol file!")
        return 1
    elif len(args) == 1:
        infile = args[0]
    else:
        parser.print_help()
        print("Too many arguments!")
        return 1
    core(infile, options)
    return 0

if __name__ == "__main__":
    #import cProfile, os
    #i = 1
    #while os.path.exists("profile_%i" % i): i += 1
    #cProfile.run('sys.exit(main(sys.argv[1:]))', "profile_%i" % i)
    sys.exit(main(sys.argv[1:]))

    # Vestigal code for validating automatic atom typing:
    #m = read_mdl_molfile(sys.argv[1])
    #add_fields_to_atoms(m.atoms)
    #assign_rosetta_types(m.atoms)
    #for i, a in enumerate(m.atoms):
    #    err_flag = ""
    #    if a.name.strip() != a.ros_type.strip():
    #        if a.name == "CH1" and a.ros_type.startswith("CH"): pass
    #        else: err_flag = "********" #raise ValueError("typing mismatch!")
    #    print "%3i %4s --> %4s %s" % (i+1, a.name, a.ros_type, err_flag)
#}}}