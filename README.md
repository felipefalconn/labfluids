# labfluids

Codes used to obtain data and perform technical analysis for experiments carried out at the Fluids Laboratory of the Federal University of Minas Gerais (UFMG), in 2025.

<img width="1360" height="655" alt="Figure_2" src="https://github.com/user-attachments/assets/764e3697-b94e-4d03-bb03-3ea985aef828" />

---

## Features

- Data processing for laboratory experiments in Fluid Mechanics  
- Computation of manometric head, uncertainties and performance curves  
- Plotting of \(H \times Q\) curves (head vs. flow rate) for different pump configurations  
- Scripts written in Python, organized per experiment

The repository is aimed at:
- Supporting lab reports
- Ensuring reproducible calculations
- Providing a clean reference implementation for future semesters

---

## Repository structure

The code is organized by experiment. A typical structure is:

```text
labfluids/
├── exp1/          # Experiment 1 (e.g. pumps in series/parallel)
│   ├── main.py    # Main script to run the analysis
│   └── ...        # Auxiliary modules, data files, figures, etc.
├── requirements.txt
└── README.md
