# OpenFF's PTM Parameterization Prototype!

## ✨ New version! ✨

The PTM prototype was updated on November 6, 2025 to be easier and more representative of our eventual goals. Changes include:

1. A new version of Pablo makes residue specification much easier
2. Simplified parametrization with the release of the OpenFF 3.0.0-alpha0 "Rosemary" Alpha force field

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

* Constructs an [openff-pablo]([url](https://github.com/openforcefield/openff-pablo/)) `ResidueDefinition` for the modified amino acid
* Solvates the modified protein in a water box
* Uses the Rosemary alpha force field to assign parameters to both canonical and non-canonical amino acids. 
* Runs a short simulation
