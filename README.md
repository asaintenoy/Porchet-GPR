# Porchet-GPR
Simulation code of radargram acquired along a Porchet water infiltration test

## Installation

- Installer miniconda

- Installer gprmax en suivant les instructions sur

  - https://github.com/gprmax/gprMax

- Installer Porchet-GPR avec

```
conda env create -f environment.yml
```

- Installer tmux

      sudo apt install tmux

- Installer dask

      sudo apt install dask

- Compiler H2D dans ./HD2

      make

## Running

```
tmux
```

Dans le nouveau shell

    python Calcul-Multi-Modeles.py

- mettre en arrière plan: `ctrl-b d`

- remettre en premier plan:

      tmux attach

D'autres infos sur tmux

- http://www.dayid.org/comp/tm.html
