# OpenFF's PTM Parameterization Prototype!

## Short description

Run this notebook to see how to parameterize proteins with PTMs using OpenFF. Then try simulating your own modified protein, and write in to this repo's issue tracker to let us know how we can make it easier for your use cases (in particular, how we can help you define `ResidueTemplates` more easily using your inputs).

## To run

```shell
git clone https://github.com/openforcefield/ptm_prototype.git
cd ptm_prototype
micromamba create -n openff_ptm_prototype -c conda-forge -f env.yaml
micromamba run -n openff_ptm_prototype jupyter lab ptm_sim.ipynb
```

## Long description

With the alpha release of "Rosemary", our upcoming protein+small molecule force field that will handle PTMs natively, we wanted to demonstrate the _infrastructure_ to do PTM parameterization now, so it will be ready on the day the final Rosemary force field is released. This notebook shows how to parameterize and simulate a modified protein using the Rosemary alpha for both canonical and noncanonical amino acids.

The notebook also demonstrates two important upcoming features that we intend to support long-term:

* [`openff-pablo`](https://github.com/openforcefield/openff-pablo): Our new, performant PDB loader that will eventually get merged into [`Topology.from_pdb`](https://docs.openforcefield.org/projects/toolkit/en/stable/users/pdb_cookbook/index.html). Currently only handles PDB files, but will handle PDBx/mmCIF in the future.
* NAGL charges: A SMIRNOFF force field using partial charges from a NAGL model assigns partial charges so quickly that the entire non-canonical protein can be charged with a consistent charge model - no library charges in sight. 

## What's in the notebook?

The `ptm_sim.ipynb` notebook in this repo shows an example of using these new features to parameterize a modified protein. As a brief overview, it:

* Constructs an [openff-pablo]([url](https://github.com/openforcefield/openff-pablo/)) `ResidueDefinition` for the modified amino acid (this is the hard part, please write into this repo's issue tracker to let us know how we can improve this).
    * In this case, we use the `ResidueDefinition.from_molecule` method, by making an OpenFF Molecule where a cysteine is attached to a covalent ligand using a SMARTS reaction.
    * You don't have to do it exactly this way for your own custom residues - We offer [different methods](https://openff-pablo.readthedocs.io/en/latest/api/generated/openff.pablo.ResidueDefinition.html) for constructing `ResidueDefinitions` and would love your suggestions for more.
    * The atom names in this `ResidueDefinition` are required to match those used for the modified residue in the input PDB.
    * If your input PDB has CONECT records for all nonstandard residues, you don't need to use openff-pablo, and can instead use the `_additional_substructures` and `_custom_substructures` arguments to `Topology.from_pdb` (this is a bit slower since it has to do graph matching, but it doesn't require you to match up atom/residue names between the PDB and substructure)
* Solvates the modified protein in a water box
* Uses the Rosemary alpha force field to assign parameters to both canonical and non-canonical amino acids. 
* Runs a short simulation

## How can I load a modified protein?

<!--TODO: Simplify this flowchart, using only Pablo for loading-->

```mermaid
flowchart TB
    0@{ shape: diamond, label: "Do you have a structure of the modified protein with all hydrogens added?" }
    0-->|No|0A
    0-->|Yes|1
    
    0A@{ shape: rectangle, label: "All OpenFF protein loading methods require explicit hydrogens to be present. Consider using PDBFixer, pdb2pqr, or another protein preparation tool."}
    0A-->|Ok, I've added explicit hydrogens|1
    
    1@{ shape: diamond, label: "Do you have a RDKit or OpenEye Molecule of your modified protein?" }
    1-->|Yes|1A
    1-->|No|2

    1A@{ shape: rectangle, label: "<pre style='white-space: pre-wrap; text-align: left; width: 400px;'>offmol = Molecule.from_rdkit(rdmol)
    # or
    offmol = Molecule.from_openeye(oemol)
    # and finally
    top = offmol.to_topology()</pre>"}
    1A-->20

    2@{ shape: diamond, label: "Do you have an SDF file of the modified protein?"}
    2-->|Yes|2A
    2-->|No|2P5
    
    2A@{ shape: diamond, label: "Do you ALSO have a PDB file with residue information that you want to preserve in the final OpenFF Molecule?"}
    2A-->|Yes|2AA
    2A-->|No|2AB


    2AA@{ shape: rectangle, label: "<pre style='white-space: pre-wrap; text-align: left; width: 1000px;'>ref_mol = Molecule.from_file('protein.sdf')
    from openmm.app import PDBFile
    pdb = PDBFile('protein.pdb')
    top = Topology.from_openmm(pdb.topology, unique_molecules=[ref_mol], positions=pdb.positions)</pre>"}
    2AA-->20
    
    2AB@{ shape: rectangle, label: "<pre style='white-space: pre-wrap; text-align: left; width: 500px;'>offmol = Molecule.from_file('protein.sdf')
    #(then optionally)
    offmol.perceive_residues()
    # and finally
    top = offmol.to_topology()</pre>"}
    2AB-->20
    
    2P5@{ shape: rectangle, label: "You must have SDF/SMILES of at least modification, if not the whole residue (for CCD residues, you can download the 'ideal' SDF from RCSB PDB)"}
    2P5-->|Yes|3
    
    3@{ shape: diamond, label: "Do you have a PDB file with CONECT records for the modified residue?"}
    3-->|Yes|4
    3-->|No|14
    
    4@{ shape: diamond, label: "Is the modified residue a conjugated cysteine?"}
    4-->|Yes|5
    4-->|No|6

    5@{ shape: rectangle, label: "<pre style='white-space: pre-wrap; text-align: left; width: 900px;'>mol = Molecule.from_file('maleimide_hydrogenated.sdf')
# Mark leaving hydrogen where the bond will be made to the cysteine
for atom in mol.atoms:
&nbsp;&nbsp;&nbsp;&nbsp;atom.metadata['substructure_atom'] = True
smarts = 'C(=O)C([H:1])CC(=O)'
# Since this pattern will be found multiple times due to symmetry, only use 
# the first match ([0]) to identify a single leaving hydrogen. 
for atom_idx in mol.chemical_environment_matches(smarts)[0]:
&nbsp;&nbsp;&nbsp;&nbsp;print(f'marking {atom_idx=} as leaving')
&nbsp;&nbsp;&nbsp;&nbsp;mol.atom(atom_idx).metadata['substructure_atom'] = False
top = Topology.from_pdb('3ip9_dye.pdb', _additional_substructures=[mol])</pre>"}
    5-->20
    
    6@{ shape: diamond, label: "Does the ENTIRE modified residue have CONECT records (not just the covalent attachment)?"}
    6-->|Yes|7
    6-->|No|14

    7@{ shape: diamond, label: "Do you have an SDF/SMILES of the ENTIRE modified residue (including backbone atoms)?"}
    7-->|Yes|10
    7-->|No|8

    8@{ shape: diamond, label: "Do you have an SDF/SMILES of the covalent attachment?"}
    8-->|Yes|9
    8-->|No|17

    9@{ shape: rectangle, label: "<pre style='white-space: pre-wrap; text-align: left; width: 1000px;'>from ptm_prototype import react
from openff.pablo import ccd, CCD_RESIDUE_DEFINITION_CACHE
# Load canonical form of residue from CCD cache
cysteine = CCD_RESIDUE_DEFINITION_CACHE['CYS'][0].to_openff_molecule()
# Load covalent attachment
maleimide = Molecule.from_file('maleimide.sdf')
# Define the reaction that conjugates the canonical residue to the covalent attachment
thiol_maleimide_click_smarts = (
    '[C:10]-[S:1]-[H:2]'
    + '.'
    + '[N:3]1-[C:4](=[O:5])-[C:6](-[H:11])=[C:7](-[H:12])-[C:8](=[O:9])-1'
    + '>>'
    + '[N:3]1-[C:4](=[O:5])-[C:6](-[H:2])(-[H:11])-[C@:7](-[S:1]-[C:10])(-[H:12])-[C:8](=[O:9])-1'
)
# Run the reaction, generate a conformer, and save the file
products = react([cysteine, maleimide], thiol_maleimide_click_smarts)
mod_res = [*products][0][0]
mod_res.generate_conformers()
mod_res.to_file('CYS-MAL.sdf', file_format='sdf')</pre>"}

    9-->10
    
    10@{ shape: rectangle, label: "<pre style='white-space: pre-wrap; text-align: left; width: 900px;'>mod_res = Molecule.from_file('CYS-MAL.sdf')
# Mark leaving hydrogen where the bond will be made to the cysteine
for atom in mod_res.atoms:
&nbsp;&nbsp;&nbsp;&nbsp;atom.metadata['substructure_atom'] = True
backbone_smarts = '[H:1][NH2X3][CH1X4][CH0X3](=[OH0X1])[OX2H1:2][H:3]'
# Since this pattern will be found multiple times due to symmetry, only use 
#the first match ([0]) to identify a single leaving hydrogen. 
for atom_idx in mod_res.chemical_environment_matches(backbone_smarts)[0]:
&nbsp;&nbsp;&nbsp;&nbsp;print(f'marking {atom_idx=} as leaving')
&nbsp;&nbsp;&nbsp;&nbsp;mod_res.atom(atom_idx).metadata['substructure_atom'] = False
top = Topology.from_pdb('3ip9_dye.pdb', _additional_substructures=[mod_res])
</pre>"}
    10-->20

    14@{ shape: rectangle, label: "Follow the instructions in our <a href='https://github.com/openforcefield/ptm_prototype' target='_blank'>PTM Prototype notebook</a>"}
    

    17@{ shape: rectangle, label: "Please write in to our <a href='https://github.com/openforcefield/ptm_prototype/issues' target='_blank'>issue tracker</a> and let us know which representation of your modified amino acid/covalent attachment you DO have."}
    17-->14

    20@{ shape: rectangle, label: "To solvate and assign parameters to the loaded Topology, start at the 'Solvate' or 'Parameterize' step in our <a href='https://github.com/openforcefield/ptm_prototype' target='_blank'>PTM Prototype notebook</a>"}


```
