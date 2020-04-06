import csv
import pandas


transcript_list = []
protein1_list = []
protein2_list = []
protein3_list = []

for i in range(1, 101):
    if i < 71:
        transcript_list.append([i, "protein1", (0.2 * i)])
    if i == 70:
        protein1_list = transcript_list[-1]
    transcript_list.append([i, "protein2", (0.12 * i)])
    if i == 100:
        protein2_list = transcript_list[-1]
    transcript_list.append([i, "protein3", (0.06 * i)])
    if i == 100:
        protein3_list = transcript_list[-1]
for i in range(71, 101):
    transcript_list.append([i, "protein1", protein1_list[-1]])
for i in range(101, 251):
    transcript_list.append([i, "protein1", protein1_list[-1]])
    transcript_list.append([i, "protein2", (0.12 * i)])
    transcript_list.append([i, "protein3", protein3_list[-1]])
# for i in range(151, 251):
#     transcript_list.append([i, "proteinX", (0.1 * i)])
#     transcript_list.append([i, "proteinY", proteinY_list[-1]])
#     transcript_list.append([i, "proteinZ", (-0.05 * i) + (0.15 * 200)])



file = open("../../data/grant_data17.tsv", "w")
writer = csv.writer(file, delimiter='\t')
writer.writerow(["time", "species", "transcript"])
for item in transcript_list:
    writer.writerow(item)
file.close()

print("Success!")
