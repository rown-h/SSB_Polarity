"""
May 2026

DESCRIPTION:
Script to export the average distance between each amino acid backbone and the 5' or 3' DNA end.

Contains one function:  dna_distances(out_dir = '', protein = '')
Optional additional arguments defined below.

USAGE:
To utilise, first prepare 'prot' and 'DNA' objects with multiple states in PyMOL

    If necessary, join multiple objects into states in one object:
        join_states states_combined, mod1 mod2 mod3 mod4 mod5...

    Separate protein to its own object with multiple states:
        create prot, states_combined and polymer.protein

    Separate DNA to its own object with multiple states:
        create DNA, states_combined and polymer.nucleic

Then run the commands in PyMOL

    Load this script:
        run /path/to/dna_distances.py

    Run the function:
        dna_distances(out_dir = '', protein = '')

ARGUMENTS:
out_dir
    string: Directory to output the CSV. Note, spaces in directory names may cause issues. {no default}

protein
    string: Name of your protein, which will be included in the filename. {default = 'my_protein'}

raw_csv
    bool: Output a csv with the individual distances between each DNA end and residue in every state. {default = False}

prot_obj
    string: Name of the isolated protein object in the PyMOL session. {default = 'prot'}

dna_obj
    string: Name of the isolated DNA object in the PyMOL session. {default = 'DNA'}

five_prime_resi
    int/string: PyMOL residue number for the 5' DNA end. If 'first', it will automatically determine this. {default = 'first'}

three_prime_resi
    int/string: PyMOL residue number for the 3' DNA end. If 'last', it will automatically determine this. {default = 'last'}

five_prime_atom
    string: PyMOL atom label from which to calculate 5' distances. If your model does not contain a complete DNA end, you may need to adjust this to "O5'". {default = "P"}

three_prime_atom
    string: PyMOL atom label from which to calculate 3' distances. Note, double quotes required if prime in atom label. {default = "O3'"}

IF ERRORS ARE OCCURING:
Check the following:
1. Check the output for the number of states detected. Are the number of states for DNA and prot the same?
2. Check the output for the number of atoms defined as the 5' phosphate selection and 3' OH selection. It should say 1 atom each. 
   If 0, set DNA to view licorice, and check whether an atom is missing. Adjust five_prime_atom if required.
3. Is the DNA oriented in the same direction in each model?
4. Manually inspect data before averaging with raw_csv = True 

"""

from pymol import cmd
import math
import csv
import os

def dna_distances(
        out_dir=None,
        protein='my_protein',
        raw_csv=False,
        prot_obj='prot',
        dna_obj='DNA',
        five_prime_resi='first',
        three_prime_resi='last',
        five_prime_atom='P',
        three_prime_atom="O3'"
        ):
    

    # Provide warnings if important arguments are not provided

    if out_dir is None:
        print("😪 Error: You need to define the full output directory path. \n          Please set out_dir = '.../SSB_Polarity/Data/Distance_Data'")
        return
    
    if protein == 'my_protein':
        print(f"🐸 Warning: Please include an argument for the name of your protein! Defaulting to protein = '{protein}'")
    
    # Concatenate to yield the out_file 
    out_file = os.path.join(out_dir, protein + '_distance.csv')

    n_states = cmd.count_states(dna_obj)
    print(f"Detected {n_states} states")

    # --- Identify 5′ and 3′ nucleotides ---
    if five_prime_resi == 'first':
        five_prime_resi = cmd.get_model(f"{dna_obj} and polymer.nucleic").atom[0].resi
    if three_prime_resi == 'last':
        model = cmd.get_model(f"{dna_obj} and polymer.nucleic")
        three_prime_resi = model.atom[-1].resi

    # Atom selections
    sel_5 = f"{dna_obj} and resi {five_prime_resi} and name {five_prime_atom}"
    sel_3 = f"{dna_obj} and resi {three_prime_resi} and name {three_prime_atom}"

    count_sel_5 = cmd.count_atoms(sel_5)
    count_sel_3 = cmd.count_atoms(sel_3)

    print(f"5′ phosphate selection: {sel_5} - {count_sel_5} atom")
    print(f"3′ OH selection: {sel_3} - {count_sel_3} atom")

    # --- Get protein residues ---
    model = cmd.get_model(f"{prot_obj} and name N")

    # Outputs tuple to identify residue, of chain, residue number, and residue name: (A, 67, LYS)
    residues = sorted(set((atom.chain, atom.resi, atom.resn) for atom in model.atom),
                      key=lambda x: int(x[1]))

    # Storage
    dist_5p = {res: [] for res in residues}
    dist_3oh = {res: [] for res in residues}

    # --- Loop over states ---
    for state in range(1, n_states + 1):

        for (chain, resi, resn) in residues:

            n_sel = f"{prot_obj} and chain {chain} and resi {resi} and name N"

            # Ensures the state is the same in both DNA and prot
            try:
                d1 = cmd.get_distance(sel_5, n_sel, state, state)
                d2 = cmd.get_distance(sel_3, n_sel, state, state)

                dist_5p[(chain, resi, resn)].append(d1)
                dist_3oh[(chain, resi, resn)].append(d2)
            
            except Exception as e:
                print(f"Error at state {state}, residue {resi}: {e}")

    # --- Statistics ---
    # Calculate mean
    def mean(vals):
        return sum(vals) / len(vals) if vals else float('nan')

    # Calculate sample standard deviation
    def sd(vals):
        if len(vals) < 2:
            return 0.0
        m = mean(vals)
        return math.sqrt(sum((x - m) ** 2 for x in vals) / (len(vals) - 1))

    # --- Write CSV ---
    with open(out_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            "Residue",
            "mean_angstrom_5",
            "sd_5",
            "mean_angstrom_3",
            "sd_3"
        ])

        for (chain, resi, resn) in residues:

            label = f"{resn}{resi}"

            writer.writerow([
                label,
                mean(dist_5p[(chain, resi, resn)]),
                sd(dist_5p[(chain, resi, resn)]),
                mean(dist_3oh[(chain, resi, resn)]),
                sd(dist_3oh[(chain, resi, resn)])
            ])

    print(f"\n📝 CSV written to: {out_file}")



    # --- Write raw data to CSV for debugging ---
    if raw_csv == True:
        raw_file = out_file.replace('.csv', '_raw.csv')

        with open(raw_file, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # header: residue + one column per state
            max_states = max(len(v) for v in dist_5p.values())
            
            header = ["Residue"]
            for i in range(1, max_states + 1):
                header.append(f"dist5_state{i}")
            for i in range(1, max_states + 1):
                header.append(f"dist3_state{i}")
            
            writer.writerow(header)

            for (chain, resi, resn) in residues:
                label = f"{resn}{resi}"
                
                vals_5 = dist_5p[(chain, resi, resn)]
                vals_3 = dist_3oh[(chain, resi, resn)]
                
                # pad shorter lists with empty values
                vals_5 = vals_5 + [''] * (max_states - len(vals_5))
                vals_3 = vals_3 + [''] * (max_states - len(vals_3))
                
                writer.writerow([label] + vals_5 + vals_3)

        print(f"📝 Raw CSV written to: {raw_file}")