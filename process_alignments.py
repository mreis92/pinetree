import sys
import csv

def process_dataset(filename):
	dataset_file = open(filename, "r")
	dataset_ids = {}
	count = 0
	
	for line in dataset_file:
		line = line.replace("\n", "")

		if line[0] == ">":
			id_seq = line[1:]
			dataset_ids[id_seq] = count
			count += 1

	dataset_file.close()
	return dataset_ids


def main():
	target_ids = process_dataset(sys.argv[1])
	miRNA_ids = process_dataset(sys.argv[2])
	
	results_file = open(sys.argv[3], 'r')
	lines = results_file.readlines()
	i = 0
	count = 0
	output = []
	
	
	while(i < len(lines)):
		line = lines[i].replace('\n', '')
		if(">>>" in line):
			miRNA = miRNA_ids[line.split(">>>")[1].split(" ")[0]]
		
		if(line.startswith(">>")):
			target_id = target_ids[line.split(" ")[0][2:]]
			
			i+=2
			
			start = int(lines[i].replace('\n', '').split("-")[2].split(":")[1])
			
			i+=3
			
			line = lines[i].replace('\n', '')
			
			for c in range(7,len(line)):
				if(line[c] != ' '):
					offset=c
					while(c < len(line) and line[c] != ' '):
						c+=1
					break
			
			miRNA_len = c-offset
			miRNA_seq = line[offset:offset+miRNA_len]
			
			i+=1
			line = lines[i].replace('\n', '')
			for c in range(offset, len(line)):
				if(line[c] != ' '):
					break
				start -= 1
				
			start = 0 if start < 0 else start
			
			i+=1
			
			line = lines[i].replace('\n', '')
			
			target_seq = line[offset:offset+miRNA_len]
			
			output.append(str(miRNA) + "," + str(target_id) + "," + str(start) + "," + miRNA_seq + "," + target_seq)
			count+=1
		
		i+=1	
		
	print count
	print "miRNA,target_id,start,miRNA_seq,target_seq"
	for line in output:
		print line

	
if __name__ == "__main__":
	main()
