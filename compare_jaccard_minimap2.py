razers_file = "test.sam"

dictionary = {}
with open(razers_file, "r") as file:
    for line in file:
        if line.startswith("@"):
            continue  # Skip headers
        fields = line.split("\t")
        read_id = fields[0]
        read_razers = read_id
        flag = int(fields[1])
        # length=int(fields[2])
        # ref_name = fields[2]
        position = int(fields[3])
        # strand = '-' if (flag & 16) else '+'
        if flag == 0:
            dictionary["@" + read_id] = position

count = 0

jaccard_file = "output.sam"
with open(jaccard_file, "r") as file:
    next(file)
    next(file)
    for line in file:
        fields = line.split("\t")
        read_id = fields[0] #.split("\t")
        read_id = read_id + ".1"
        read_jaccard = read_id
        flag = int(fields[1])
        length = int(fields[2])
        position = int(fields[4])
        if read_id in dictionary:
            pos = dictionary[read_id]
            if position >= (pos - length) and position <= (pos + length):
                count = count + 1

percent = (count / len(dictionary)) * 100
print(f"Accuracy: {percent:.2f}%")
print(count)

print(read_razers)
print(read_jaccard)
