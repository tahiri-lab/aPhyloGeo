reference: data
		python reference.py

data: fetch_data.sh
		./fetch_data.sh

.PHONY: clean

clean: 
		rm output/*