data=[]
f = open("data.txt", 'r') 
line=f.readline()
while line:
    result=line.split()
    data.append(result)
    line=f.readline()
print len(data)
 
Adata = []
with open('Adata.txt') as f:
 
    for line in f:                   # loop over the rows
        fields = line.split()        # parse the columns
        rowdata = map(float, fields) # convert text to numbers
        Adata.extend(rowdata)         # accumulate the results
# print 'A Minimum:', min(Adata), 'A Maximum:', max(Adata)
 
Cdata = []
with open('Cdata.txt') as f:
    for line in f:                   # loop over the rows
        fields = line.split()        # parse the columns
        rowdata = map(float, fields) # convert text to numbers
        Cdata.extend(rowdata)         # accumulate the results
# print len(data), "                                                        #length of data set"
print min(Adata), max(Adata), min(Cdata), max(Cdata), " #x_min x_max y_min y_max"
 
#print min(Adata), max(Adata), min(Cdata), max(Cdata), " #x_min x_max y_min y_max"
  
edata = []
with open('edata.txt') as f:
   for line in f:                   # loop over the rows
       fields = line.split()        # parse the columns
       rowdata = map(float, fields) # convert text to numbers
       edata.extend(rowdata)         # accumulate the results
# print len(data), "
#edata = edata[:] - min(data)
# openf.write( len(data) )
# openf.write(min(Adata), max(Adata), min(Cdata), max(Cdata), " #x_min x_max y_min y_max")
# print data
#nmine=min(edata)
for line in range(len(data)):
    print data[line][0], data[line][1], float(data[line][2])-float(min(edata))*0.999999
#   print data[line][0], data[line][1], data[line][2]
# openf.write(data[line][0], data[line][1], data[line][2])
