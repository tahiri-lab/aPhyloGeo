aPhylogeo:
	python aphylogeo/main.py run

.PHONY: clean

clean:
	rm -rf output/*
	rm -rf *_newick
	rm -rf output.csv

test:
	python -m pip install -r requirements-dev.txt
	python -m pytest tests
