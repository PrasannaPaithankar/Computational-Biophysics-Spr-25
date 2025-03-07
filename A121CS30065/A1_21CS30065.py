# Assignment 1: Protein Sequence and Structure Analysis
# Course: Computational Biophysics: Algorithms to Applications (CS61060)
# Author: Prasanna Paithankar
# Roll Number: 21CS30065
# Date: 07-03-2025

import os
import sys
from datetime import datetime

import requests
import toml
from Bio import PDB, SeqIO, Align


def output_header():
    return f"""------------------------------------------------------------------------
Assignment 1: Protein Sequence and Structure Analysis
Course: Computational Biophysics: Algorithms to Applications (CS61060)
Author: Prasanna Paithankar
Roll Number: 21CS30065
Generated on: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
------------------------------------------------------------------------
"""


def download_file(ID: str, file_type: str) -> str:
    """
    Downloads PDB or FASTA file for a given ID
    Stores the file in the cache directory
    """

    if file_type == "pdb":
        url = f"{config['PDB_BASE_URL']}{ID}.pdb"
    elif file_type == "fasta":
        url = f"{config['FASTA_BASE_URL']}{ID}/download"
    else:
        raise Exception(f"Unsupported file type: {file_type}")

    response = requests.get(url)
    if response.status_code == 200:
        with open(f"cache/{ID}.{file_type}", "wb") as file:
            file.write(response.content)
        return f"cache/{ID}.{file_type}"
    else:
        raise Exception(f"Failed to download {file_type} file for {ID}")


def extract_fasta_sequence(fasta_file: str) -> dict:
    sequences = {}
    with open(fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            chains_t = record.description.split("|")[1]
            if chains_t.startswith("Chains"):
                chains = chains_t[7:]
                chains = chains.split(", ")
            else:
                chains = chains_t[6:]
            for chain in chains:
                sequences[chain] = str(record.seq)
    sequences = dict(sorted(sequences.items()))
    return sequences


def extract_pdb_sequence(pdb_file: str) -> dict:
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("PDB", pdb_file)
    sequences = {}
    for model in structure:
        for chain in model:
            chain_id = chain.id
            sequence = ""
            for residue in chain:
                if residue.get_id()[0] == " ":
                    resname = residue.get_resname()
                    one_letter_code = config["AMINO_ACID_CODES"][resname]
                    if one_letter_code:
                        sequence += one_letter_code
            sequences[chain_id] = sequence
    return sequences


def identify_chain_breaks(pdb_chains: dict, fasta_sequence: dict) -> dict:
    aligner = Align.PairwiseAligner()
    chain_breaks = {}

    for chain_id, pdb_seq in pdb_chains.items():
        fasta_seq = fasta_sequence.get(chain_id, "")
        alignments = aligner.align(pdb_seq, fasta_seq)[0]
        rel2fasta = alignments.aligned[1]
        rel2fasta = [(int(a), int(b)) for (a, b) in rel2fasta]
        breaks = []
        start = 0
        for (a, b) in rel2fasta:
            if a > start:
                breaks.append(f"({start + 1} to {a})")
            start = b
        if start < len(fasta_seq):
            breaks.append(f"({start + 1} to {len(fasta_seq)})")

        if breaks:
            chain_breaks[chain_id] = breaks
        else:
            chain_breaks[chain_id] = None

    return chain_breaks


def naccess(ID: str, pdb_file: str) -> None:
    """
    The NACCESS path can be set in the config.toml file
    """

    os.system(f"{config['NACCESS_PATH']} {pdb_file}")
    os.rename(f"{ID}.asa", f"cache/{ID}.asa")
    os.rename(f"{ID}.rsa", f"cache/{ID}.rsa")
    os.rename(f"{ID}.log", f"cache/{ID}.log")


def parse_asa(ID: str, pdb_file: str) -> dict:
    naccess(ID, pdb_file)
    asa = {}
    with open(f"cache/{ID}.rsa", "r") as file:
        for line in file:
            if line.startswith("CHAIN"):
                chain_id = line.split()[2]
                asa[chain_id] = float(line.split()[3])
    return asa


def compute_molecular_weight(sequence: str) -> float:
    weights = config["AMINO_ACID_WEIGHTS"]
    return sum(weights[aa] for aa in sequence)


def main(ID: str) -> None:
    try:
        if not os.path.exists("cache"):
            os.makedirs("cache")
        pdb_file = f"cache/{ID}.pdb"
        fasta_file = f"cache/{ID}.fasta"

        # Download files if not present in cache
        if not (os.path.exists(pdb_file) and os.path.exists(fasta_file)):
            pdb_file = download_file(ID, "pdb")
            fasta_file = download_file(ID, "fasta")

        # Extract sequences
        fasta_sequence = extract_fasta_sequence(fasta_file)
        pdb_chains = extract_pdb_sequence(pdb_file)

        # Identify chain breaks
        chain_breaks = identify_chain_breaks(pdb_chains, fasta_sequence)

        # Parse the NACCESS output to get accessible surface area (ASA)
        asa = parse_asa(ID, pdb_file)

        # Create report file
        with open(f"A1_{ID}.txt", "w") as file:
            file.write(output_header())
            file.write(f"\nID: {ID}\n")
            file.write(f"\nNumber of chains: {len(pdb_chains)}\n")
            file.write("\nFASTA sequences:\n")
            for chain_id, sequence in fasta_sequence.items():
                file.write(f"Chain {chain_id}: {sequence}\n")
            file.write("\nPDB chains:\n")
            for chain_id, sequence in pdb_chains.items():
                file.write(f"Chain {chain_id}: {sequence}\n")
            file.write("\nChain breaks:\tFormat = [(break1_start to break1_end), (break2_start to brak2_end), ...]\tBoth ends are inclusive\n")
            for chain_id, breaks in chain_breaks.items():
                    if breaks:
                        file.write(f"Chain {chain_id}: {breaks}\n")
                    else:
                        file.write(f"Chain {chain_id}: No chain breaks\n")
            file.write(
                "\n{:<10} {:<20} {:<20} {:<10}\n".format(
                    "Chain ID", "No. of Amino Acids", "Molecular Weight", "ASA"
                )
            )
            for chain_id, sequence in pdb_chains.items():
                file.write(
                    "{:<10} {:<20} {:<20.2f} {:<10.2f}\n".format(
                        chain_id,
                        len(sequence),
                        compute_molecular_weight(sequence),
                        asa.get(chain_id, 0),
                    )
                )

    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    ID = "1UBQ"
    config = toml.load("config.toml")
    if len(sys.argv) > 1:
        ID = sys.argv[1]
    else:
        ID = config["ID"]

    main(ID)
