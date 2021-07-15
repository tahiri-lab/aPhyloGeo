reference: data
		python3 reference.py

matrix: input_matrix.txt
		python3 matrix.py

data: fetch_data.sh
		./fetch_data.sh

.PHONY: clean

clean: 
		rm output/*