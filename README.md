# Welcome to the OpenFF PTM Parameterization Prototype!

Since "Rosemary" (our upcoming protein+small molecule force field that will handle PTMs natively) is taking a while, we wanted to start the ball rolling on the infrastructure to do PTM parameterization NOW, by using a chimera of the SMIRNOFF port of AMBER's FF14SB and OpenFF's Sage. You may have seen experimental demos of this sort of workflow in our previous virtual workshops, but in contrast to those one-and-done notebooks, this prototype demonstrates the use of two upcoming features that we intend to support long-term:

* `openff-pablo`: Our new, performant polymer loader that will eventually supersede `Topology.from_pdb`. Currently only handles PDB files, but will handle PDBx/mmCIF in the future.
* `parameterize_with_nagl`: A concept of the "swiss cheese" partial charge assignment method, where one force field (in this case, our FF14SB port) does as much as it can to assign partial charges to a molecule, and then a NAGL model  assigns the rest. 

The `ptm_sim.ipynb` notebook in this repo shows our recommended use of these new features to parameterize a modified protein. As a brief overview, it:

* Constructs a new residue definition for the modified amino acid by conjugating canonical cysteine to a small molecule using a SMARTS reaction.
    * Requires some work to line these atom names up with those in the PDB residue (this is the hard part)
    * If the element+connectivity based template matching already available in Topology.from_pdb is sufficient for you (namely, if your input PDB file has CONECT records) you can skip this
* Solvates the modified protein in a water box
* Assigns canonical residues ff14SB parameters, and everything else gets Sage parameters with NAGL charges
* Runs a short simulation


To run:

```shell
git clone https://github.com/openforcefield/ptm_prototype.git
cd ptm_prototype
micromamba create -n openff_ptm_prototype -c conda-forge openff-toolkit-examples python=3.11 pyxdg
micromamba run -n openff_ptm_prototype pip install git+https://github.com/openforcefield/openff-pablo.git@v0.0.1a1
micromamba run -n openff_ptm_prototype jupyter lab ptm_sim.ipynb
```

Distributions like Arch with very recent Glibc (>=2.41) don't play well with the version of PyTorch installed in the above configuration. If you get an "ImportError: libtorch_cpu.so: cannot enable executable stack as shared object requires: Invalid argument" error, try:

```shell
git clone https://github.com/openforcefield/ptm_prototype.git
cd ptm_prototype
micromamba create -n openff_ptm_prototype -f env.yaml
micromamba run -n openff_ptm_prototype jupyter lab ptm_sim.ipynb
```
