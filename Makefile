reference: scripts/fetch_data.sh scripts/names.py pipeline.py
		./scripts/fetch_data.sh
		python3 scripts/names.py
		python3 pipeline.py

matrix: input_files/input.txt
		python3 tree.py

.PHONY: clean

clean: 
		rm output/windows/*