import csv
import os
for i in os.listdir("."):
	if i.endswith(".tsv"):
		with open(i) as f, open(i+".len","w") as o:
    			reader = csv.reader(f, delimiter='\t')
    			writer = csv.writer(o, delimiter='\t')
    			for row in reader:
        			writer.writerow([row[0],len(row[1])])
