import csv
import pandas


transcript_list = []
protein1_list = []
protein2_list = []
protein3_list = []

for i in range(1, 251):
    transcript_list.append([i, "protein1", (0.2 * i)])
    transcript_list.append([i, "protein2", (0.12 * i)])
    transcript_list.append([i, "protein3", (0.04 * i)])
    # transcript_list.append([i, "protein4", (0.14 * i)])
    # transcript_list.append([i, "protein5", (0.12 * i)])
    # transcript_list.append([i, "protein6", (0.10 * i)])
    # transcript_list.append([i, "protein7", (0.08 * i)])
    # transcript_list.append([i, "protein8", (0.06 * i)])
    # transcript_list.append([i, "protein9", (0.04 * i)])
    # transcript_list.append([i, "protein10", (0.02 * i)])
# for i in range(71, 101):
#     transcript_list.append([i, "protein1", protein1_list[-1]])
# for i in range(101, 251):
#     transcript_list.append([i, "protein1", protein1_list[-1]])
#     transcript_list.append([i, "protein2", (0.12 * i)])
#     transcript_list.append([i, "protein3", protein3_list[-1]])
# for i in range(151, 251):
#     transcript_list.append([i, "proteinX", (0.1 * i)])
#     transcript_list.append([i, "proteinY", proteinY_list[-1]])
#     transcript_list.append([i, "proteinZ", (-0.05 * i) + (0.15 * 200)])



file = open("../../data/grant_data19.tsv", "w")
writer = csv.writer(file, delimiter='\t')
writer.writerow(["time", "species", "transcript"])
for item in transcript_list:
    writer.writerow(item)
file.close()

print("Success!")
