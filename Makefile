reference: data
		python3 pipeline.py

matrix: input_files/input.txt
		python3 matrix.py

data: fetch_data.sh
		./fetch_data.sh

.PHONY: clean

clean: 
		rm output/windows/*