import json
import pandas as pd
import matplotlib.pyplot as plt


with open("./daily.json", "r") as json_file:
    daily = json.load(json_file)

Vx = []
Vy = []
Vz = []
Vs = [Vx, Vy, Vz]
for n, (T, meter) in enumerate(daily.items()):
    V = Vs[n]
    for m, data in meter.items():
        for day, cols in data.items():
            V_data = cols["Voltage"]
            V += V_data


plt.plot(Vs[0])
plt.plot(Vs[1])
plt.plot(Vs[2])
plt.title("Perfiles Datos")
plt.ylabel("|V| [V]")
plt.xlabel("Muestra")
plt.legend(("TRA21089",
           "TRA1624",
           "TRA24198"))
plt.savefig("./Vprofile.pdf")
plt.show()