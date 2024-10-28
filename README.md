[![CI](https://github.com/delphes/delphes/actions/workflows/ci.yml/badge.svg)](https://github.com/delphes/delphes/actions/workflows/ci.yml) [![DOI](https://zenodo.org/badge/21390046.svg)](https://zenodo.org/badge/latestdoi/21390046)

Delphes
=======

Delphes is a C++ framework, performing a fast multipurpose detector response simulation.

More details can be found on the Delphes website http://cp3.irmp.ucl.ac.be/projects/delphes

Quick start with Delphes
========================

Commands to get the code:

```
  git clone git@github.com:MuonColliderSoft/Delphes.git Delphes

  cd Delphes

  source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc12-opt/setup.sh

  make
```

The card for the MUSIC detector can be found in `Delphes/cards/delphes_card_MUSICDet_target.tcl`

To perform some validation of the MUSIC Delphes card, you would need first to build Delphes with Pythia:

Build Pythia
```
wget https://pythia.org/download/pythia83/pythia8310.tgz
tar xzvf pythia8310.tgz
cd pythia8310
./configure --prefix=path_to_PYTHIA8_installation
make install
```

Second, define an environment variable for the path to your PYTHIA installation directory
```
export PYTHIA8=path_to_PYTHIA8_installation
```
and you can then build the DelphesPythia8 executable with the following command:
```
make HAS_PYTHIA8=true
```

Then, you can go inside the `validation` folder and follow the instructions in the README file.