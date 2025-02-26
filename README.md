# Competitive Binding Simulator
Python based simulator for competitive binding experiments, e.g. fluorescence polarization. 
Equations taken from: Michael H. A. Roehrl, Julia Y. Wang, and Gerhard Wagner, Biochemistry 43, 16056 (2004). 

Depends on: matplotlib, numpy

```
usage: competitive_binding_sim.py [-h] [-o OUTPUT] [-kd KD] [-ki KI]
                                  [-rt RECEPTOR_TOTAL] [-lt LIGAND_TOTAL]

Plot simulated competitive binding curves for the input values

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output filename for plot
  -kd KD, --kd KD       Binding constant for the labeled ligand (micromolar)
  -ki KI, --ki KI       Binding constant for the competitor (micromolar)
  -rt RECEPTOR_TOTAL, --receptor_total RECEPTOR_TOTAL
                        Concentration of receptor (normally protein;
                        micromolar)
  -lt LIGAND_TOTAL, --ligand_total LIGAND_TOTAL
                        Concentration of labeled ligand; micromolar
```

Notes: 

 - It is assumed that competition is "complete"; i.e. binding of the labeled ligand and competitor are mutually exclusive.
 - No sanity checks are performed on inputs at present.
 - Units for all input parameters are micromolar.
