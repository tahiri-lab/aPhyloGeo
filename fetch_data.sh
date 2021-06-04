#!/bin/bash

# This script retrieves all the fasta files in the data repository and create a new fasta file named 'sequences.fasta' containing all the sequences.
# This file will then serves as an input in the muscle executable.

# Set the filename
filename='sequences.fasta'

# Check if the file already exists, if it does, delete it
if [ -f $filename ]; then
    rm sequences.fasta
fi

# Create a new sequences.fasta file containing all the sequences in the data repo
find data -name "*.fasta" -exec cat {} >> output/sequences.fasta \;

