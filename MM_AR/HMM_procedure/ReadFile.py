import csv

def ReadDataPoint(datafile):
    '''Generator object creating GPS Point'''
    csvfile = open(datafile, 'rb')
   
    reader = csv.reader(csvfile)              
    rownum = 0       
    for row in reader:
        # Save header row.
        if rownum == 0:
            #header = row
            pass
        else:
            #colnum = 0
            #for col in row:
            #order in csv file 
            #for 20140307.csv     ['lat', 'lon', 'ele', 'time', 'course', 'speed', 'src', 'sat', 'hdop']
            
            yield (row[5], row[4],row[9])
            
            #colnum += 1                
        rownum += 1
        
    csvfile.close()
    
    #yield (gps_lon, gps_lat, timestamp)

#datafile = r'C:\Users\asr13006\Google Drive\UConn MS\Data Collection\test20140218_upto154700.csv'
#datafile = datafile.replace("\\", "/")
#print datafile

#for i in ReadDataPoint(datafile):
#    print i