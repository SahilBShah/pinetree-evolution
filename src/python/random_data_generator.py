import csv
import pandas


transcript_list = []
proteinX_list = []
proteinY_list = []
proteinZ_list = []

for i in range(1, 81):
    transcript_list.append([i, "proteinX", (0.3 * i)])
    if i == 80:
        proteinX_list = transcript_list[-1]
    transcript_list.append([i, "proteinY", (0.4 * i)])
    if i == 80:
        proteinY_list = transcript_list[-1]
    transcript_list.append([i, "proteinZ", (0.2 * i)])
    if i == 80:
        proteinZ_list = transcript_list[-1]
for i in range(81, 151):
    transcript_list.append([i, "proteinX", proteinX_list[-1]])
    transcript_list.append([i, "proteinY", proteinY_list[-1]])
    transcript_list.append([i, "proteinZ", proteinZ_list[-1]])
for i in range(151, 241):
    transcript_list.append([i, "proteinX", proteinX_list[-1]])
    transcript_list.append([i, "proteinY", (-0.15 * i) + (0.4 * 80) + (0.15 * 151)])
    transcript_list.append([i, "proteinZ", proteinZ_list[-1]])



file = open("../../data/grant_data/grant_data10.tsv", "w")
writer = csv.writer(file, delimiter='\t')
writer.writerow(["time", "species", "transcript"])
for item in transcript_list:
    writer.writerow(item)
file.close()

print("Success!")
