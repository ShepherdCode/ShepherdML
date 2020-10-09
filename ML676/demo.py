import csv
FILENAME = "Cars.csv"
mycars = []
with open (FILENAME) as handle:
    reader = csv.DictReader(handle)
    for line in reader:
        mycars.append(line)
for car in mycars:
    if int(car['PriceRangeLow'])<=30000 and int(car['PriceRangeHigh'])>=30000:
        print(car['Model'],car['PriceRangeLow'])
