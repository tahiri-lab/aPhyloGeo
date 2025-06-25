import pandas as pd
import numpy as np
from sklearn.feature_selection import VarianceThreshold
from Bio import AlignIO
from io import StringIO
import json

def filter_low_variance_features(csv_path, threshold=0.01, id_column='id'):
    """
    Filter the columns with low variance in a CSV file.

    Parameters:
    - csv_path: Route to the CSV file.
    - threshold: Minimum variance threshold.
    - id_column: Column that identifies the specimens.

    Returns:
    - filtered_df: DataFrame with only the useful columns.
    """
    df = pd.read_csv(csv_path)
    
    # Separe the ID column from the rest of the DataFrame
    if id_column in df.columns:
        metadata = df[[id_column]]
        features = df.drop(columns=[id_column])
    else:
        metadata = pd.DataFrame()
        features = df.copy()
    
    # Aplies the variance filter
    selector = VarianceThreshold(threshold)
    filtered_array = selector.fit_transform(features)
    selected_columns = features.columns[selector.get_support(indices=True)]

    # Rebuild the DataFrame with useful columns
    filtered_df = pd.concat([metadata.reset_index(drop=True),
                             pd.DataFrame(filtered_array, columns=selected_columns)], axis=1)
    return filtered_df

def preprocess_windowed_alignment(input_path, threshold=0.8, output_path = None ):
    """
    Remove the columns with too many gaps in each window of the alignment, preserving the window format.

    Parameters:
    - input_path: file JSON with 'msa': {window: alignment in FASTA format}
    - output_path: resulting JSON file
    - threshold: maximum ratio of gaps allowed per column (default 0.8)
    """
    with open(input_path, "r") as f:
        data = json.load(f)

    original_msa = data.get("msa", {})
    filtered_msa = {}

    for window, fasta_str in original_msa.items():
        alignment = AlignIO.read(StringIO(fasta_str), "fasta")
        keep_columns = []

        for i in range(alignment.get_alignment_length()):
            column = alignment[:, i]
            gap_count = column.count("-") + column.upper().count("N")
            if (gap_count / len(column)) <= threshold:
                keep_columns.append(i)

        # We build new sequences with only the allowed columns
        new_seqs = {}
        for record in alignment:
            new_seq = "".join([record.seq[i] for i in keep_columns])
            new_seqs[record.id] = new_seq

        # We rebuild the FASTA in string
        filtered_fasta = "\n".join(f">{seq_id}\n{seq}" for seq_id, seq in new_seqs.items())
        filtered_msa[window] = filtered_fasta

    # We save a new file with the same format but filtered windows
    with open(output_path, "w") as f:
        json.dump({"type": "Alignment", "alignment_method": "1", "msa": filtered_msa}, f, indent=2)


# === Previous versions for the preprocessing ====
"""
def preprocess_alignment_remove_columns_with_gaps(file_path, threshold=0.1, output_path=None):
    with open(file_path) as f:
        data = json.load(f)

    # Detecting if the format is {id: seq} or {"msa": {...}}
    if "msa" in data:
        full_sequences = {}
        for window in data["msa"].values():
            for entry in window.strip().split(">"):
                if entry.strip() == "":
                    continue
                lines = entry.strip().split("\n")
                seq_id = lines[0]
                seq = "".join(lines[1:])
                full_sequences[seq_id] = full_sequences.get(seq_id, "") + seq
    else:
        # Direct format {id: seq}
        full_sequences = data

    # Transform the matrix (list of lists of characters)
    ids = list(full_sequences.keys())
    sequences = [list(full_sequences[_id]) for _id in ids]
    sequence_length = len(sequences[0])
    num_sequences = len(sequences)

    # Verifiying that all sequences have the same length
    for seq in sequences:
        if len(seq) != sequence_length:
            raise ValueError("Las secuencias no tienen la misma longitud")

    # Detecting columns with too many gaps
    columns_to_keep = [
        col for col in range(sequence_length)
        if sum(seq[col] == '-' for seq in sequences) / num_sequences <= threshold
    ]

    # Build new sequences with only the columns to keep
    filtered_sequences = {
        _id: "".join([sequences[i][j] for j in columns_to_keep])
        for i, _id in enumerate(ids)
    }

    if output_path:
        with open(output_path, "w") as f:
            json.dump(filtered_sequences, f, indent=2)

    return filtered_sequences
"""

"""
def filter_genetic_windows_by_gaps(window, threshold=0.2):
    
    Filters out a window of sequences if any sequence exceeds the allowed proportion of gaps ('-' or 'N').

    Args:
        window (dict): Dictionary of sequences {id: sequence}.
        threshold (float): Maximum allowed proportion of gaps per sequence.

    Returns:
        dict: Same structure as `window`, only if all sequences pass the filter.
              Returns an empty dictionary if any sequence exceeds the threshold.
    
    for seq_id, seq in window.items():
        gap_count = seq.count('-') + seq.upper().count('N')
        gap_ratio = gap_count / len(seq)
        print(f"{seq_id}: {gap_ratio:.2%} gaps")
        if gap_ratio > threshold:
            print(f"Sequence {seq_id} exceeds the gap threshold ({threshold:.0%})")
            return {}  # Discard the entire window
    print("Window accepted")
    return window



def filter_low_variance_matrix(matrix, threshold=0.01):

    Filters rows and columns with low variance from a symmetric matrix (e.g., genetic distance matrix).

    Parameters:
    - matrix: DataFrame or numpy array (square matrix)
    - threshold: Minimum variance threshold

    Returns:
    - filtered_matrix: The filtered symmetric matrix
    - kept_indices: Indices of rows/columns retained

    if isinstance(matrix, pd.DataFrame):
        data = matrix.values
    else:
        data = matrix

    # Calculate the variance for each column (which will be the same for each row if it is symmetric)
    variances = np.var(data, axis=0)
    kept_indices = np.where(variances >= threshold)[0]

    #Filter the matrix keeping the symmetry
    filtered_matrix = data[np.ix_(kept_indices, kept_indices)]

    #If it was a DataFrame, return as DataFrame
    if isinstance(matrix, pd.DataFrame):
        kept_labels = matrix.columns[kept_indices]
        return pd.DataFrame(filtered_matrix, index=kept_labels, columns=kept_labels), kept_indices
    else:
        return filtered_matrix, kept_indices
"""
