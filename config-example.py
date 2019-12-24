import datetime

host='192.168.0.69'
database='aurora-test'
user='root'
password='root'

# threshold boundaries evaluating, relative to maximum flux power value
threshold = 0.1

# set datetimeToProcess to None for next forecast datetime auto evaluating
datetimeToProcess = datetime.datetime(2019,4,1, 0,0,0) # 1 april 2019, 00:00