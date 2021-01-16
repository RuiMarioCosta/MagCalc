# MagCalc
MagCalc calculates finite temperatures properties such as: magnetization, entropies (magnetic, lattice and total), free energies and others as functions of temperature and applied magnetic field.

## requirements
* python3
* numpy
* scipy
* matplotlib
* pyyaml

## usage
By executing the program with:

`python main.py`

a prompted will appear to choose which of the materials defined in `magcalc/configuration.yaml` you want to use.

Alternatively the name of the material can be passed as an argument to `--name` which skips the prompt, e.g.,

`python main.py --name Gd5Si2Ge2`