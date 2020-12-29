import os

for filename in os.listdir("genomes"):
    f = open("genomes/"+filename, "r")
    total_len = 0
    for line in f.readlines()[1:]:
        total_len += len(line)
    print("file name: {:33}  file length: {}".format(filename, total_len))
    f.close()

