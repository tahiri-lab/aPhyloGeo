import re
import os

def test():
    f = open("infile", "r")
    line_count = -1
    for line in f:
        if line != "\n":
            line_count += 1
    f.close()
    f = open("infile", "r").read()
    num_seq = int((f.split("\n")[0]).split(" ")[0]) # get the number of sequences in this file
    longueur = int((f.split("\n")[0]).split(" ")[1])
    no_line = int(line_count/num_seq)
    with open("outfile", "w") as out:
        depart = 1
        fin = depart + no_line 
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
    
    # open outfile and align all the sequences on top of another
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
    with open("out", "r") as f, open("outfile", "w") as out:
        for line in f:
            out.write(line[0:20] + "\n") # HERE, we should change the value corresponding to the steps

        f.close()
        out.close()
    with open("out", "w") as f, open("outfile", "r") as out:
        f.write(str(num_seq) + " " + str(longueur) + "\n")
        for index, line in enumerate(out):
            if line != "\n":
                f.write(list_names[index] + "\t")
                f.write(line)

        


if __name__ == "__main__":
    test()
