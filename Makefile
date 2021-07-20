reference: fetch_data.sh pipeline.py
		./fetch_data.sh
		python3 pipeline.py

matrix: input_files/input.txt
		python3 tree.py

.PHONY: clean

clean: 
		rm output/windows/*