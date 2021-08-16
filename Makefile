reference: scripts/fetch_data.sh scripts/names.py pipeline.py
		./scripts/fetch_data.sh
		python3 scripts/names.py
		python3 pipeline.py < in.txt > test

tree: input/input.txt
		python3 tree.py

.PHONY: clean

clean:
		rm output/reference_gene.fasta