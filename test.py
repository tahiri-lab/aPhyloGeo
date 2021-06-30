import re
import os

def test(window_size=0, step=0):
    # Permet d'avoir le nombre de lignes totales dans le fichier
    f = open("infile", "r")
    line_count = -1
    for line in f:
        if line != "\n":
            line_count += 1
    f.close()
    f = open("infile", "r").read()
    num_seq = int((f.split("\n")[0]).split(" ")[0]) # premier nombre de la premiere ligne du fichier represente le nbr de sequences
    longueur = int((f.split("\n")[0]).split(" ")[1]) # second nombre de la premiere ligne du fichier represente la longueur des sequences
    no_line = int(line_count/num_seq) # permet d'obtenir le nbr de lignes qui compose chaque sequence

    # Recupere la sequence pour chaque variante
    with open("outfile", "w") as out:
        depart = 1
        fin = depart + no_line 
        # on connait la longueur de chaque sequence, donc on va recuperer chaque sequence et le retranscrire sur un autre fichier separes par un \n entre chaque
        for i in range(0, int(num_seq)):
            f = open("infile", "r")
            lines_to_read = range(depart,fin)
            for position,line in enumerate(f):
                if position in lines_to_read:
                    out.write(line)
            out.write("\n")
            depart = fin
            fin = depart + no_line
    out.close()
    
    # on cree un fichier out qui contient chaque sequence sans espaces et on enregistre dans une list le nom en ordre des sequences
    with open("outfile", "r") as out, open("out", "w") as f:
        sequences = out.read().split("\n\n")
        list_names = []
        for seq in sequences:
            s = seq.replace("\n", " ").split(" ")
            if s[0] != "":
                list_names.append(s[0])
            s_line = s[1:len(seq)]
            for line in s_line:
                if line != "":
                    f.write(line)
            f.write("\n")
    out.close()
    f.close()
    
    # slide the window along the sequence
    debut = 0
    fin = debut + window_size
    while fin <= longueur:
        index = 0
        with open("out", "r") as f, open("output/windows/" + str(debut), "w") as out:
            out.write(str(num_seq) + " " + str(window_size) + "\n")
            for line in f:
                if line != "\n":
                    espece = list_names[index]
                    nbr_espaces = 11 - len(espece)
                    out.write(espece)
                    for i in range(nbr_espaces):
                        out.write(" ")
                    out.write(line[debut:fin] + "\n") 
                    index = index + 1
        out.close()
        f.close()
        debut = debut + step
        fin = fin + step

        


if __name__ == "__main__":
    test(10,3)
