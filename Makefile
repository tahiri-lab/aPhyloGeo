aPhylogeo: 
	python3 ./scripts/main.py

.PHONY: clean

clean:
	rm -rf output/*
	rm -rf *_newick
	rm -rf output.csv
