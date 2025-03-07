### 🧬 Computational Biophysics: Algorithms to Applications (CS61060)
# Assignment 1: Protein Sequence and Structure Analysis
#### 🧑🏻‍🔬 Prasanna Paithankar (21CS30065)

## 📝 Problem Objectives
1. Download a PDB structure and the fasta sequence file from RCSB Protein Data Bank (rcsb.org).
2. Extract the protein sequence from the PDB file and compare it with the fasta sequence file.
3. Identify the number of chains in the protein structure.
4. Check if there are any chain breaks in the structure or not .
5. Run the NACCESS program to compute the solvent accessible surface area (ASA) of chains.

## 🧪 Usage
```bash
uv run A1_21CS30065.py <PDB ID>
```

#### The **NACCESS path** is stored in `config.toml` file.

**Note**: The PDB ID can also be stored in `config.toml` file.

## 📤 Output
A output file named `A1_{PDB ID}.txt` will be generated in the same directory.

***
#### 2025 | Department of Computer Science and Engineering, IIT Kharagpur
