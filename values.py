from calculations import *
import csv

data={
    "experimental":"$4464.980\\pm0.20$",
    "ground-state-energy":ground_state_energy(),
    "expectation-of-the-delta":expectation_of_the_delta(),
    "hyperfine-splitting":hyperfine_splitting()
    }

with open('../Computing-project-report/data.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["name", "value"])
    for key, value in data.items():
        writer.writerow([key, value])