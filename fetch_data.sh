#!/bin/bash

# This script retrieves all the fasta files in the data repository and create a new fasta file named 'sequences.fasta' containing all the sequences.
# This file will then serves as an input in the muscle executable.

find data -name *.fasta -exec cat {} >> sequences.fasta \;

