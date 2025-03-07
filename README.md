# OpenFF's PTM Parameterization Prototype!

## Short description

Run this notebook to see how to parameterize proteins with PTMs using OpenFF. Then try simulating your own modified protein, and write in to this repo's ssue tracker to let us know how we can make it easier for your use cases (in particular, how we can help you define `ResidueTemplates` more easily using your inputs).

## To run

```shell
git clone https://github.com/openforcefield/ptm_prototype.git
cd ptm_prototype
micromamba create -n openff_ptm_prototype -c conda-forge openff-toolkit-examples python=3.11 pyxdg
micromamba run -n openff_ptm_prototype pip install git+https://github.com/openforcefield/openff-pablo.git@v0.0.1a1
micromamba run -n openff_ptm_prototype jupyter lab ptm_sim.ipynb
```

## Long description

Since "Rosemary" (our upcoming protein+small molecule force field that will handle PTMs natively) is taking a while, we wanted to start building the _infrastructure_ to do PTM parameterization now, so it will be ready on the day that Rosemary force field is released. This notebook shows how to parameterize and simulate a modified protein using a chimera of the SMIRNOFF port of AMBER's FF14SB for canonical amino acids and OpenFF's Sage for noncanonical amino acids.

The notebook also demonstrates two important upcoming features that we intend to support long-term:

* `openff-pablo`: Our new, performant PDB loader that will eventually supersede `Topology.from_pdb`. Currently only handles PDB files, but will handle PDBx/mmCIF in the future.
* `parameterize_with_nagl`: A concept of the "swiss cheese" partial charge assignment method, where one force field with library charges (in this case, our FF14SB port) assigns partial charges to as much of a molecule as it can, and then a NAGL model assigns partial charges to the rest. 

Once Rosemary comes out, the `parameterize_with_nagl` function will likely become unnecessary (since Rosemary is on track to use a NAGL model for all charge assignment), and we will switch to recommending Rosemary alone instead of FF14SB+Sage. However the loading functionality in `openff-pablo` will continue to be needed to load modified proteins into the OpenFF ecosystem.

## What's in the notebook?

The `ptm_sim.ipynb` notebook in this repo shows an example of using these new features to parameterize a modified protein. As a brief overview, it:

* Constructs an [openff-pablo]([url](https://github.com/openforcefield/openff-pablo/)) `ResidueDefinition` for the modified amino acid (this is the hard part, please write into this repo's issue tracker to let us know how we can improve this).
    * In this case, we use the `ResidueDefinition.from_molecule` method, by making an OpenFF Molecule where a cysteine is attached to a covalent ligand using a SMARTS reaction.
    * You don't have to do it exactly this way for your own custom residues - We offer [different methods](https://openff-pablo.readthedocs.io/en/latest/api/generated/openff.pablo.ResidueDefinition.html) for constructing `ResidueDefinitions` and would love your suggestions for more.
    * The atom names in this `ResidueDefinition` are required to match those used for the modified residue in the input PDB.
    * If the element+connectivity based template matching that's already available in `Topology.from_pdb` is sufficient for you (namely, if your input PDB file has CONECT records) you can continue using that, and you don't need to make an openff-pablo `ResidueDefinition`.
* Solvates the modified protein in a water box
* Uses the `parameterize_with_nagl` method to assign ff14SB parameters to unmodified amino acids, and Sage parameters with NAGL charges to everything else. 
* Runs a short simulation

## Debugging occasional environment issues

Distributions like Arch with very recent Glibc (>=2.41) don't play well with the version of PyTorch installed in the above configuration. If you get an "ImportError: libtorch_cpu.so: cannot enable executable stack as shared object requires: Invalid argument" error, try:

```shell
git clone https://github.com/openforcefield/ptm_prototype.git
cd ptm_prototype
micromamba create -n openff_ptm_prototype -f env.yaml
micromamba run -n openff_ptm_prototype jupyter lab ptm_sim.ipynb
```
