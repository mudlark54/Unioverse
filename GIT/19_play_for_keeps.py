#------------------------------------
# read the whole file into a list
# read EVERY OTHER line, if it has the _R, then print out the NEXT Line

fout1 = open("output1.txt", 'w') #open file for writin
with open("Pfile.fas") as fin1:
	lines = fin1.readlines() #reads all of the lines into a list []
 	for i in range(0, len(lines), 2): #standard for loop for iterating over the INDEX (start at 0, go to lenght of list, count by 2)
		if lines[i].strip().endswith(("_R", "__comp0")): #lines[i] gets the line at that index. if it has the string you are looking for
											#NOTE: .strip() removes the witespace (includeing end-of-line characters). 
											#without this, the string actually ends with '\n' and the .endswith("_R") would fail
			fout1.writelines(lines[i])
			fout1.writelines(lines[i+1])


#------------------------------------
#read the file line by line. 
#if you find the _R - read the next line and print it out
#if you dont, read the next line and dump it

fin2 = open("Pfile.fas")  #open file for reading
fout2 = open("output2.txt", 'w') #open file for writing
line = fin2.readline() #reads the FIRST LINE - the first name
while line: #this will continue unil the next readline() returns the end of file
	gene = fin2.readline() #read the next line, will be a gene
	if line.strip().endswith(("_R", "__comp0")):  #endswith can take a tuple (a, b, c, d) and will return true if the string ends with ANY of them
		fout2.writelines(line)
		fout2.writelines(gene)
	line = fin2.readline() # read the next line, will be a name


#------------------------------------
# read the whole file into a list


fin3 = open("Pfile.fas") #open file for reading
fout3 = open("output3.txt", 'w') #open file for writin
lines = fin3.readlines() #reads all of the lines into a list []
names = lines[::2] #grabs every other item (starting with the first element [index=0])
genes = lines[1::2] #grabs every other ittem (**starrting with the second elementt [index = 1])
for i in range(len(names)):
	if names[i].strip().endswith(("_R", "__comp0")):
		fout3.writelines(names[i])
		fout3.writelines(genes[i])


#------------------------------------


fin4 = open("Pfile.fas") #open file for reading
fout4 = open("output4.txt", 'w') #open file for writin
lines = fin4.readlines() #reads all of the lines into a list []
names = lines[::2] 
genes = lines[1::2] 
nameGeneDict = zip(names, genes) #creates a dictinary by merging 2 lists {name:gene, name:gene}
for x,y in nameGeneDict:
	if x.strip().endswith(("_R", "__comp0")):
		fout4.writelines(x)
		fout4.writelines(y)


#------------------------------------


#List comprehention (my favorite) 1 liner!
with open("output5.txt", 'w') as fout: fout.writelines([x+y for x, y in zip(open("Pfile.fas").readlines()[::2], open("Pfile.fas").readlines()[1::2]) if x.strip().endswith(("_R", "__comp0"))])


#------------------------------------

