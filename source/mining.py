import numpy as np
from Bio.PDB import PDBParser, Superimposer, Selection, PDBIO, Select
from Bio.SeqUtils import seq1
from pathlib import Path

# if you want to use the debugger, uncomment the following line
# from IPython.core.debugger import set_trace

from scipy.spatial import cKDTree
import subprocess
import os
import argparse
import sys


def load_structure(file_path):
    parser = PDBParser(QUIET=True)
    return parser.get_structure(file_path.stem, file_path)


def parse_dssp(dssp_fn):
    secondary_structure = {}
    with open(dssp_fn, "r") as f:
        for line in f:
            if line.startswith("  #  RESIDUE"):
                break
        for line in f:
            if line[13] != "!":
                res_id = int(line[5:10])
                ss = line[16]
                secondary_structure[res_id] = ss
    return secondary_structure


def run_dssp(structure_fn, stem_fn, glycine_res_id, temp_output_dir):
    dssp_executable = "mkdssp"
    dssp_fn = temp_output_dir + stem_fn.stem + "_matched_{}.dssp".format(glycine_res_id)
    command = [dssp_executable, "-i", structure_fn, "-o", dssp_fn]

    dssp_result = subprocess.run(command)

    secondary_structure_dssp = parse_dssp(dssp_fn)

    secondary_structure = [
        secondary_structure_dssp.get(i, "?")
        for i in range(glycine_res_id - 6, glycine_res_id + 5)
    ]
    os.remove(dssp_fn)
    return secondary_structure


def find_motifs(residues, motif_length, glycine_pos):
    """Yields continuous motifs of a given length that have a glycine at a specific position."""

    for i in range(len(residues) - motif_length + 1):
        segment = residues[i : i + motif_length]
        if segment[glycine_pos].get_resname() == "GLY":
            yield segment


def calc_rmsd(seg1, seg2):
    """Calculates the RMSD between two segments based on backbone atoms."""
    sup = Superimposer()
    atoms1 = (
        [r["N"] for r in seg1 if "N" in r.child_dict]
        + [r["CA"] for r in seg1 if "CA" in r.child_dict]
        + [r["C"] for r in seg1 if "C" in r.child_dict]
    )
    atoms2 = (
        [r["N"] for r in seg2 if "N" in r.child_dict]
        + [r["CA"] for r in seg2 if "CA" in r.child_dict]
        + [r["C"] for r in seg2 if "C" in r.child_dict]
    )

    if len(atoms1) == len(atoms2):
        sup.set_atoms(atoms2, atoms1)
        return sup.rms, sup
    else:
        return float("inf"), None


def align_structure(sup, structure):
    """Applies the superimposer transformation to a segment"""
    # Apply transformation to each atom in the motif
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom.transform(sup.rotran[0], sup.rotran[1])


def count_clashes_to_crbn(structure, crbn_structure_pcd_tree, radius=2.0):
    """Counts the number of heavy atom and ca clashes between the structure and the CRBN structure."""

    structure_coords = [
        atom.get_coord() for atom in structure.get_atoms() if atom.element != "H"
    ]
    structure_ca_coords = [
        atom.get_coord() for atom in structure.get_atoms() if atom.get_id() == "CA"
    ]

    # Invoke cKDTree with (crbn_structure_pcd_tree)
    d_nn_ca, i_nn_ca = crbn_structure_pcd_tree.query(
        structure_ca_coords, k=1, distance_upper_bound=radius
    )
    d_nn, i_nn = crbn_structure_pcd_tree.query(
        structure_coords, k=1, distance_upper_bound=radius
    )
    clashing_ca = np.sum(d_nn_ca <= radius)
    clashing = np.sum(d_nn <= radius)
    return clashing_ca, clashing


class NoHydrogensNoDisorder(Select):
    def accept_atom(self, atom):
        if atom.element == "H":
            return False
        else:
            if (not atom.is_disordered()) or atom.get_altloc() == "A":
                atom.set_altloc(" ")  # Eliminate alt location ID before output.
                return True
            else:
                return False


def save_structure(matched_structure, file_path):
    io = PDBIO()
    with open(file_path, "w") as f:
        io.set_structure(matched_structure)
        io.save(f, NoHydrogensNoDisorder())


def add_cryst1_line(file_path):
    # Add the CRYST1 line at the beginning.
    line = "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1"
    with open(file_path, "r") as original:
        data = original.read()
    with open(file_path, "w") as modified:
        modified.write(line + "\n" + data)


def main(
    query_fn_list,
    crbn_structure_path,
    glycine_pos,
    template_motif_fn,
    output_dir,
    temp_output_dir,
    pdb_dir=None,
    af2_dir=None,
    rmsd_threshold=1.0,
):

    print("*************")
    print("Query file list:", query_fn_list)

    template_motif_struct = load_structure(Path(template_motif_fn))
    template_motif = Selection.unfold_entities(template_motif_struct, "R")

    crbn_structure = load_structure(Path(crbn_structure_path))
    crbn_structure_coords = [
        atom.get_coord() for atom in crbn_structure.get_atoms() if atom.element != "H"
    ]
    crbn_structure_pcd_tree = cKDTree(crbn_structure_coords)

    out_csv_fn = (
        output_dir + "/" + os.path.basename(query_fn_list).replace(".list", ".csv")
    )

    out_csv = open(out_csv_fn, "w")
    template_motif_length = len(template_motif)
    out_csv.write("PDBID,RMSD,Clashing_CA,Clashing,")
    for i in range(template_motif_length):
        out_csv.write(f"Res{i},")
    out_csv.write("Degron,SS\n")

    with open(query_fn_list, "r") as f:
        for line in f.readlines():
            fn = line.replace("pdb_domains", pdb_dir).strip()
            fn = fn.replace("af2_domains", af2_dir).strip()
            fn = Path(fn)
            print("Processing:", fn)
            if not fn.is_file() or fn.suffix != ".pdb":
                continue
            print("Matching motifs in:", fn.stem)
            query_structure = load_structure(fn)

            # Remove hydrogens from structure
            query_structure_residues = Selection.unfold_entities(query_structure, "R")

            for matched_motif in find_motifs(
                query_structure_residues, len(template_motif), glycine_pos
            ):

                rmsd, sup = calc_rmsd(matched_motif, template_motif)

                if rmsd <= rmsd_threshold:
                    structure_copy = query_structure.copy()

                    try:
                        for model in structure_copy:
                            for chain in model:
                                chain.id = "B"

                    except Exception:
                        print("%% Error changing chain id for %s" % fn.stem)
                        continue

                    # Align the matched motif (and the rest of the structure) to the template motif
                    align_structure(sup, structure_copy)
                    clashing_ca, clashing = count_clashes_to_crbn(
                        structure_copy, crbn_structure_pcd_tree
                    )
                    if clashing_ca <= 1 and clashing <= 12:
                        print(f"Motif found with RMSD: {rmsd} for {fn.stem}")
                        pdbid = fn.stem

                        print(f"Clashing CA: {clashing_ca}, Clashing: {clashing}")

                        outline = f"{pdbid},{rmsd},{clashing_ca},{clashing}"

                        glycine_res_id = matched_motif[glycine_pos].get_id()[1]

                        for res in matched_motif:
                            outline += f",{res.get_resname()} {res.get_id()[1]}"
                            print(f"{res.get_resname()} {res.get_id()}", end=" ")

                        outdegron = ""
                        # Get the 14 residues before and the 14 residues after the glycine
                        for i in range(glycine_res_id - 14, glycine_res_id + 14):
                            if (" ", i, " ") in structure_copy[0]["B"]:
                                resname = structure_copy[0]["B"][i].get_resname()
                                oneletter = seq1(resname)
                            else:
                                oneletter = "X"
                            outdegron += oneletter
                        outline += f",{outdegron}"

                        # Save as a new PDB file
                        new_fn = (
                            temp_output_dir + fn.stem + f"_matched_{glycine_res_id}.pdb"
                        )
                        save_structure(structure_copy, new_fn)

                        # concatenate new_fn and crbn_structure_path into a new file
                        filenames = [crbn_structure_path, new_fn]
                        concat_fn = temp_output_dir + fn.stem + "_concat.pdb"
                        with open(concat_fn, "w") as outfile:
                            for fname in filenames:
                                with open(fname) as infile:
                                    for line in infile:
                                        outfile.write(line)

                        secondary_structure = ""
                        # Run DSSP on the matched file after adding the cryst1 line
                        try:
                            add_cryst1_line(new_fn)
                            secondary_structure = run_dssp(new_fn, fn, glycine_res_id, temp_output_dir)
                        except Exception:
                            print("%%Error running DSSP for %s" % fn.stem)
                            # print the actual error and exception
                            print(sys.exc_info())

                            continue

                        outline += f",{''.join(secondary_structure)}\n"
                        out_csv.write(outline)

                        # Remove the temporary files
                        os.remove(concat_fn)
                        # os.remove(new_fn) # this is the aligned structure, keep it for testing purposes

                    print("\n########")
            print("###############################")


if __name__ == "__main__":
    # Example usage
    crbn_structure_path = "input_template/CRBNx5FQDxB.pdb"
    # Use argparse to parse arguments, take as argument query_file_list, glycine_pos, target_motif_pdb fn

    parser = argparse.ArgumentParser(description="Process mining parameters..")

    # Add arguments
    parser.add_argument(
        "query_file_list",
        type=str,
        help="Path to the file containing a list of query files.",
    )
    parser.add_argument(
        "--glycine_pos",
        type=int,
        default=5,
        help="Position of glycine in the sequence.",
    )
    parser.add_argument(
        "--pdb_dir",
        type=str,
        help="Path to the .pdb files.",
    )
    parser.add_argument(
        "--af2_dir",
        type=str,
        default=None,
        help="Path to the AlphaFold2 files.",
    )
    parser.add_argument(
        "--target_motif_pdb_fn",
        type=str,
        default=None,
        help="Filename of the target motif in PDB format.",
    )
    parser.add_argument(
        "--output_dir", type=str, required=True, help="Directory to store output files."
    )

    parser.add_argument(
        "--temp_output_dir",
        type=str,
        required=True, 
        help="Directory to store temporary alignment files."
    )

    # Parse the arguments
    args = parser.parse_args()

    # Here you could add functionality to use these arguments as needed

    main(
        query_fn_list= args.query_file_list,
        crbn_structure_path = crbn_structure_path,
        glycine_pos = args.glycine_pos,
        template_motif_fn= args.target_motif_pdb_fn,
        pdb_dir = args.pdb_dir,
        af2_dir = args.af2_dir,
        output_dir = args.output_dir,
        temp_output_dir = args.temp_output_dir
    )
