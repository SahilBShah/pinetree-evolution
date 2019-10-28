import csv
import pandas


transcript_list = []
proteinX_list = []
proteinY_list = []
proteinZ_list = []

for i in range(1, 151):
    transcript_list.append([i, "proteinX", (0.1 * i)])
    if i == 150:
        proteinX_list = transcript_list[-1]
    transcript_list.append([i, "proteinY", (0.2 * i)])
    if i == 150:
        proteinY_list = transcript_list[-1]
    transcript_list.append([i, "proteinZ", (0.15 * i)])
    if i == 150:
        proteinZ_list = transcript_list[-1]
for i in range(151, 251):
    transcript_list.append([i, "proteinX", (0.1 * i)])
    transcript_list.append([i, "proteinY", proteinY_list[-1]])
    transcript_list.append([i, "proteinZ", (-0.05 * i) + (0.15 * 200)])



file = open("../../data/grant_data/grant_data13.tsv", "w")
writer = csv.writer(file, delimiter='\t')
writer.writerow(["time", "species", "transcript"])
for item in transcript_list:
    writer.writerow(item)
file.close()

print("Success!")
