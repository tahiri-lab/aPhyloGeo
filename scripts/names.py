sequences_file = open("output/reference_gene.fasta", "r")
list_of_lines = sequences_file.readlines()
for index in range(len(list_of_lines)):
    if list_of_lines[index].startswith(">"):
        splitted_line = list_of_lines[index].split("/")
        name = ">" + splitted_line[2] + "\n"
        list_of_lines[index] = name

sequences_file = open("output/reference_gene.fasta", "w")
sequences_file.writelines(list_of_lines)
sequences_file.close()
