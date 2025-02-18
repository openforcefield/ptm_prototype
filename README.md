# ptm_prototype

To run:

```shell
git clone https://github.com/openforcefield/ptm_prototype.git
cd ptm_prototype
micromamba create -n openff_ptm_prototype openff-toolkit-examples python=3.11 pyxdg
micromamba run -n openff_ptm_prototype pip install git+https://github.com/openforcefield/openff-pablo.git@v0.0.1a1
micromamba run -n openff_ptm_prototype jupyter lab ptm_sim.ipynb
```

Distributions like Arch with very recent Glibc (<=2.41) don't play well with the version of PyTorch installed in the above configuration. If you get an "ImportError: libtorch_cpu.so: cannot enable executable stack as shared object requires: Invalid argument" error, try:

```shell
git clone https://github.com/openforcefield/ptm_prototype.git
cd ptm_prototype
micromamba create -n openff_ptm_prototype -f env.yaml
micromamba run -n openff_ptm_prototype jupyter lab ptm_sim.ipynb
```
