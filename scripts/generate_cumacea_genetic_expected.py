import os
from pathlib import Path
from io import StringIO
from Bio import AlignIO, Phylo
from aphylogeo import utils
from aphylogeo.alignement import AlignSequences
from aphylogeo.params import Params
import json

# === Load Params and sequences ===
Params.load_from_file("tests/cumacea_align.yaml")
ref_gene_dir = Params.reference_gene_dir
ref_gene_file = Params.reference_gene_file
sequences = utils.loadSequenceFile(os.path.join(ref_gene_dir, ref_gene_file))
seq_alignment = AlignSequences(sequences.copy())

test_dir = Path("tests/testFiles")
dataset_name = "Cumacea"

print(f"bootstrap_threshold: {Params.bootstrap_threshold}")
print(f"rate_similarity: {Params.rate_similarity}")
print(f"dist_threshold: {Params.dist_threshold}")
print(f"window_size: {Params.window_size}")
print(f"step_size: {Params.step_size}")
print(f"alignment_method: {Params.alignment_method}")
print(f"tree_type: {Params.tree_type}")

# === Centroid ===
centroid = seq_alignment.getSequenceCentroid()[0]
centroid_dir = test_dir / "getSequenceCentroid" / dataset_name
centroid_dir.parent.mkdir(parents=True, exist_ok=True)
with open(centroid_dir, "w") as f:
    f.write(str(centroid))
print(f"[Centroid] {centroid}")

# === Pairwise alignments ===
aligned = seq_alignment.alignSequencesWithPairwise(centroid, sequences.pop(centroid))
align_dir = test_dir / "alignSequence" / dataset_name
align_dir.mkdir(parents=True, exist_ok=True)
for key, aln in aligned.items():
    seq_alignment.dictToFile(aln, align_dir / f"{key}.fasta", ".fasta")
print(f"[Pairwise Alignments] n_alignments: {len(aligned)}")

# === Pad sequences to same length ===
def pad_sequences_to_max_length(seq_dict):
    lengths = [len(seq) for seq in seq_dict.values()]
    if len(set(lengths)) == 1:
        return seq_dict
    max_len = max(lengths)
    return {k: seq.ljust(max_len, "-") for k, seq in seq_dict.items()}

# === Sliding window to FASTA string blocks ===
def split_alignment_into_windows(msa_dict, window_size, step_size):
    windows = {}
    sequence_length = len(next(iter(msa_dict.values())))
    for start in range(0, sequence_length, step_size):
        end = start + window_size
        if end > sequence_length:
            break
        label = f"{start}_{end - 1}"
        fasta_block = ""
        for seq_id, seq in msa_dict.items():
            fasta_block += f">{seq_id}\n{seq[start:end]}\n"
        windows[label] = fasta_block
    return windows

# === Heuristic MSA ===
heuristicMSA = seq_alignment.starAlignement(centroid, aligned)
heuristicMSA = pad_sequences_to_max_length(heuristicMSA)

# === Save windowed alignment as JSON ===
windowed_msa = split_alignment_into_windows(
    heuristicMSA,
    window_size=Params.window_size,
    step_size=Params.step_size
)
windowed_output_path = test_dir / "windowed_alignment" / f"{dataset_name}.json"
windowed_output_path.parent.mkdir(parents=True, exist_ok=True)
utils.dictToJson(windowed_msa, windowed_output_path)
print(f"[Windowed MSA] saved to {windowed_output_path}")

# === Save full heuristic MSA ===
staralign_dir = test_dir / "starAlignement" / dataset_name
staralign_dir.mkdir(parents=True, exist_ok=True)
seq_alignment.dictToFile(heuristicMSA, staralign_dir / f"{dataset_name}.fasta", ".fasta")
print(f"[Heuristic MSA] sequences: {len(heuristicMSA)}")

# === Sliding Window MSA ===
windowed = seq_alignment.slidingWindow(heuristicMSA)
genetic_list = list(windowed.keys())
genetic_list_dir = test_dir / "createGeneticList"
genetic_list_dir.mkdir(parents=True, exist_ok=True)
with open(genetic_list_dir / f"{dataset_name}.txt", "w") as f:
    f.write(str(genetic_list))

sliding_dir = test_dir / "slidingWindow" / dataset_name
sliding_dir.mkdir(parents=True, exist_ok=True)
for key, aln in windowed.items():
    seq_alignment.dictToFile(aln, sliding_dir / f"{key}.fasta", ".fasta")
print(f"[Sliding Window] n_windows: {len(genetic_list)} | windows: {genetic_list}")

# === MSA por ventana ===
msa = seq_alignment.makeMSA(windowed)
msa_dir = test_dir / "makeMSA" / dataset_name
msa_dir.mkdir(parents=True, exist_ok=True)
for key, aln in msa.items():
    with open(msa_dir / f"{key}.fasta", "w") as f:
        AlignIO.write(aln, f, "fasta")
print(f"[MSA] n_MSA: {len(msa)}")

# === FastTree Trees ===
fasttree_dir = test_dir / "fasttree" / dataset_name
fasttree_dir.mkdir(parents=True, exist_ok=True)
if Params.tree_type == '2':
    fasttree_trees = utils.fasttree(msa, Params.bootstrap_threshold, True)
    for key, tree in fasttree_trees.items():
        with open(fasttree_dir / f"tree_{key}.xml", "w") as f:
            Phylo.write(tree, f, "phyloxml")
    print(f"[FastTree] n_trees: {len(fasttree_trees)}")

# === Bootstrap Trees (esto restaurado) ===
create_bootstrap_dir = test_dir / "createBootstrap" / dataset_name
create_bootstrap_dir.mkdir(parents=True, exist_ok=True)
genetic_trees = utils.createBoostrap(msa, Params.bootstrap_threshold)
with open(create_bootstrap_dir / f"{dataset_name}.xml", "w") as f:
    for tree in genetic_trees.values():
        Phylo.write(tree, f, "phyloxml")

# Guardar bootstrapList
bootstrap_list = []
for key, tree in genetic_trees.items():
    confidences = [clade.confidence for clade in tree.find_clades() if hasattr(clade, "confidence") and clade.confidence is not None]
    bootstrap_list.append(max(confidences) if confidences else None)

bootstrap_list_dir = test_dir / "bootstrapList"
bootstrap_list_dir.mkdir(parents=True, exist_ok=True)
with open(bootstrap_list_dir / f"{dataset_name}.txt", "w") as f:
    f.write(str(bootstrap_list))
print(f"[Bootstrap Trees] n_trees: {len(genetic_trees)}")
if bootstrap_list:
    print(f"[Bootstrap List] first value: {bootstrap_list[0]}")

# === Export Cumacea JSON-formatted windowed MSA ===
cumacea_json = {
    "type": "Alignment",
    "alignment_method": str(Params.alignment_method),
    "msa": {}
}

for window_label, alignment in msa.items():
    fasta_io = StringIO()
    AlignIO.write(alignment, fasta_io, "fasta")
    cumacea_json["msa"][window_label] = fasta_io.getvalue()

os.makedirs("results/cumacea", exist_ok=True)
with open("results/cumacea/aligned_sequences_cumacea.fasta.json", "w") as f:
    json.dump(cumacea_json, f, indent=2)
print("[✓] File 'aligned_sequences_cumacea.fasta.json' generated at: results/cumacea/")
