from calculations import *
import csv

'''
data={
    "experimental":"$4464.980\\pm0.020$",
    "frolov":"$$"
    "ground-state-energy":ground_state_energy(),
    "expectation-of-the-delta":expectation_of_the_delta(),
    "hyperfine-splitting":hyperfine_splitting()
    }'
'''

data={
    "energy_Hartrees":"$"+str(energy_Hartrees())+"$",
    "energy_KeV":"$"+str(energy_KeV()).replace("+/-","\\pm")+"$"
    }

with open('../Computing-project-report/data.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["name", "value"])
    for key, value in data.items():
        writer.writerow([key, value])
        print(key,value)