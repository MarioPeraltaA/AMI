"""kWh load consumption.

Author: Mario R. Peralta A.

"""

kWhDays = []
with open("./MAR_LoadsLV.dss", "r") as file:
    for line in file:
        if "kWh" in line:
            j = line.index("!")
            i = j - 10
            val = line[i:j]
            val = val.strip()
            if "=" in val:
                val = val.split("=")[1]
            val = (float(val)) / (30*24)
            line = line.replace("kWh=", f"kWh={val}")
            j = line.index("!") + 1
            i = j - 5
            line = line.replace(line[i:j], " !")
            kWhDays.append(line)


with open("./kWhDays.dss", "w") as file:
    for line in kWhDays:
        file.write(line)
