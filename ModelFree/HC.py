"""Data Processor.

This program retrieve, sort and process data
from the utility company
as well as manages Geographic Information System (GIS).

Author: Mario R. Peralta A.
Email: Mario.Peralta@ucr.ac.cr

"""



import glob    # To deal with directories
import pandas as pd       # To call data
import geopandas as gpd   # GIS
# To creat points and lines in HIS
from shapely.geometry import Point
import folium   # Customize map
import json     # To write data as Json


class Profile():

    def __init__(self):
        pass

    def call_loads(
            self,
            path: str = "./Perfiles/Puntos_de_medicion.xlsx"
    ) -> dict:
        """Retrieve costumers data.

        It reads the network Measurement Points
        data and set it as dictionaries attributes
        one per excel sheet as follow:

            :py:attr:`Layers._loads`
            :py:attr:`Layers._meters`

        """
        # Empty structure data new attribute
        self._loads = {}
        self._meters = {}
        # Read file
        df = pd.read_excel(path,
                           sheet_name=None)
        # List of all sheets
        sheets = list(df.keys())

        # Set data regarding the sheet
        for sheet in df.keys():
            # Measurement points
            if sheet == sheets[0]:
                for c in df[sheet].columns:
                    values = [v for v in df[sheet][c]]
                    # Update attribute
                    self._loads[c] = values
            # 741: Meters information
            elif sheet == sheets[1]:
                for c in df[sheet].columns:
                    values = [v for v in df[sheet][c]]
                    # Update attribute
                    self._meters[c] = values

        return self._loads

    def rename_meter(self) -> dict:
        # Rename meters ID
        load_class = self._loads["Tarifa"]
        old_names = self._loads["Medidor"]
        new_names = {}
        for (n, load) in enumerate(load_class):
            name = old_names[n]
            if "Residencial" in load:
                new_name = f"{name}R"
            elif "Comercio" in load:
                new_name = f"{name}C"
            elif "Industria" in load:
                new_name = f"{name}I"

            new_names[name] = new_name

        self._class = new_names        # New attribute
        return new_names

    def add_loads_layer(self) -> gpd.GeoDataFrame:
        """GIS objects.

        It creats geodataframe out of :py:class:`Layers` objects
        attribute.

        Note: Devices ["ME202487", "ME264803"] are not in layers
        :py:attr:`Layers.meters_gdf` but in file data:
    
            "./Perfiles/Puntos_de_medicion.xlsx"

        However there is no profile data to the first load with
        deviceID: "ME202487"

        Note: It sets the attribute :py:attr:`Layers.loads_gdf`.

        """
        loads = self.call_loads()
        long = loads["X"]
        lat = loads["Y"]

        geometry = [Point(x, y) for x, y in zip(long, lat)]
        # Dict to df
        df = pd.DataFrame.from_dict(loads)
        # Dataframe to gdf
        gdf = gpd.GeoDataFrame(df, geometry=geometry, crs="EPSG:5367")

        return gdf

    def set_loads_attributes(self) -> gpd.GeoDataFrame:
        loads_layer = self.add_loads_layer()
        kWhmonth = self.set_monthly_Edemand()
        objID = self.rename_meter()
        # Set loads ID based on devices
        meters = [objID[m] for m in loads_layer["Medidor"]]
        loads_layer["ICEobjID"] = meters
        loads_layer["CLUSTER"] = loads_layer["Transformador"]
        loads_layer["CLASS"] = [c[-1] for c in meters]
        loads_layer["NomVoltage"] = [float(v.split("_")[1])
                                     for v in loads_layer["Tensión"]]
        loads_layer["NOMVOLT"] = loads_layer["NomVoltage"].apply(set_NOMVOLT)
        loads_layer["SERVICE"] = loads_layer["Número de fases"].apply(set_service)
        # Drop devices whose profile data is missing
        i_row = []
        for i, c in enumerate(loads_layer["ICEobjID"]):
            if c not in kWhmonth:
                i_row.append(i)
        loads_layer.drop(i_row, inplace=True)
        loads_layer["KWHMONTH"] = [kWhmonth[c]["KWHMONTH"]
                                   for c in loads_layer["ICEobjID"]]
        # Retrieve required attributes only
        attrs = [
            "CLUSTER",
            "geometry",
            "ICEobjID",
            "CLASS",
            "NomVoltage",
            "NOMVOLT",
            "SERVICE",
            "KWHMONTH",
            "X",
            "Y"
        ]
        return loads_layer[attrs]

    def read_voltage_profile(
            self,
            path: str = "./Perfiles/Perfil_de_tension.xlsx"
    ) -> dict:
        """Retrieve voltage data.

        It calls the instantaneous voltage measured data.
        It creats attribute :py:attr:`Profile._voltage`

        """
        v_meters = {}
        v_df = pd.read_excel(path,
                             sheet_name=0)

        meters = list(set(v_df["MEDIDOR"]))
        for m in meters:
            costumer = {}
            for i in v_df["MEDIDOR"].loc[v_df["MEDIDOR"] == m].index:
                time = v_df["FECHA"][i]
                voltage = v_df["INST_VOLT_A"][i]
                costumer[time] = voltage
            v_meters[m] = costumer

        # Creat new attribute
        self._voltage = v_meters

        return v_meters

    def read_load_profile(
            self,
            path: str = "./Perfiles/loadshape.xlsx"
    ) -> dict:
        """Retrieve loadshape data.

        It calls the instantaneous power demanded
        measured data.
        It creats attribute :py:attr:`Profile._load`

        """
        c_meters = {}    # Costumer meters
        c_df = pd.read_excel(path,
                             sheet_name=0)

        meters = list(set(c_df["MEDIDOR"]))
        # Creat new colum for each instantaneous
        c_df["Time"] = c_df["FECHA"] + " " + c_df["HORA"]
        for m in meters:
            costumer = {}
            for i in c_df["MEDIDOR"].loc[c_df["MEDIDOR"] == m].index:
                time = c_df["Time"][i]
                P_kW = c_df["KW"][i]
                Q_kVAr = c_df["KVAR"][i]
                S_kVA = c_df["KVA"][i]
                # Sort as active | reactive | apparent
                power = [P_kW, Q_kVAr, S_kVA]
                costumer[time] = power
            c_meters[m] = costumer

        # Creat new attribute
        self._load = c_meters

        return c_meters

    def zip_profiles(self) -> tuple[dict]:
        """Map profiles.

        Set instantaneous data available in both profiles
        in order to stablish one on one relationship between them.
        In case they fall apart with span of 5min a interpolation
        will be assessed otherwise would not be considered;
        besides, it creates attributes :py:attr:`Profile._demand_curves`
        and :py:attr:`Profile._voltage_curves` for only those voltages
        that followed such criteria.

        Note: Costumer "ME202487" is not in load curve
        but in voltage profile connected at "TRA24198"

        Note: Creates :py:attr:`Profile._odd_costumer` attribute
        which stores the smart meter ID of those costumer
        missing in the voltage profile.

        """
        self._odd_costumer = []
        loads = self._load
        voltages = self._voltage   # Profile
        # Instantaneous in both profiles
        demand_curves = {}
        voltage_curves = {}
        for m, time in loads.items():
            power_t = {}
            voltages_t = {}
            try:
                v_date = voltages[m]
            except KeyError as e:
                self._odd_costumer.append(e)
                print(f"MissData: Meter {e} does not exist.")

            for t in time.keys():
                if t in v_date:
                    voltage = v_date[t]
                    voltages_t[t] = voltage
                    power_t[t] = loads[m][t]
                else:
                    ts = pd.to_datetime(t, format='%d-%m-%Y %H:%M:%S')
                    dt = pd.Timedelta('5m')
                    ti = ts + dt
                    tj = ts - dt
                    label_i = ti.strftime("%d-%m-%Y %H:%M:%S")
                    label_j = tj.strftime("%d-%m-%Y %H:%M:%S")
                    if (label_i in v_date) and (label_j in v_date):
                        # Get average
                        vi = v_date[label_i]
                        vj = v_date[label_j]
                        Vavg = (vi+vj) / 2
                        voltages_t[t] = Vavg
                        power_t[t] = loads[m][t]
                    else:
                        continue
            voltage_curves[m] = voltages_t
            demand_curves[m] = power_t

        self._voltage_curves = voltage_curves
        self._demand_curves = demand_curves
        return (voltage_curves, demand_curves)

    def sort_profile(self) -> dict[dict[list]]:
        """Map profiles.

        It uses the time (date) as sorted index
        in order to allocate each voltage measure to each load
        power measure in a one on one link.

        Note: It creates the attribute :py:attr:`Profile._profiles`.

        """
        voltages, loads = self.zip_profiles()
        meters_ind = {}
        profiles = {}
        # Get sorted index dates
        for m, time in loads.items():
            dates = []   # Time
            for t in time.keys():
                dates.append(t)
            dates.sort()
            meters_ind[m] = dates

        for m in loads.keys():
            # Instantaneous values
            v_t = []
            p_t = []
            q_t = []
            s_t = []
            for t in meters_ind[m]:
                try:
                    v = float(voltages[m][t])   # V
                    p = float(loads[m][t][0])   # kW
                    q = float(loads[m][t][1])   # kVAr
                    s = float(loads[m][t][2])   # kVA
                    v_t.append(v)
                    p_t.append(p)
                    q_t.append(q)
                    s_t.append(s)
                except KeyError as e:
                    print(f"MissingData: {m} {e}")
                    return False
            # Sorted profiles
            data = {
                "Time": meters_ind[m],
                "Voltage": v_t,
                "Active": p_t,
                "Reactive": q_t,
                "Apparent": s_t
            }

            profiles[m] = data

        self._profiles = profiles
        return profiles

    def write_json(self, data: dict, path: str) -> None:
        """Write Json file.

        Save json file in current directory.

        """
        with open(path, "w") as file:
            # To get utf8-encoded file instead of ascii-encoded
            json.dump(data, file, ensure_ascii=False)

    def count_daily_samples(self) -> dict[dict[int]]:
        """Number of samples each day per costumer.

        Counter of daily voltage measurements.

        """
        voltage_curves, _ = self.zip_profiles()
        # All possible days
        days = set()
        for v in voltage_curves.values():
            for t in v.keys():
                day = t.split(" ")[0]
                days.add(day)
        days = list(days)
        days.sort()

        # Initialize counter and data structure
        daily_samples = {}
        for m in voltage_curves.keys():
            samples = {}
            for day in days:
                samples[day] = 0
            daily_samples[m] = samples

        # Iterate over data
        for m, data in voltage_curves.items():
            for date in data.keys():
                day = date.split(" ")[0]
                daily_samples[m][day] += 1

        return daily_samples

    def get_iday(self,
                 min_samples: int = 85) -> dict[list[tuple]]:
        """Integer index of filtered days.

        It retrieves the indices of those days with
        at least ``min_samples`` measurements
        throughout the day.

        """
        daily_i = {}
        n_samples = self.count_daily_samples()
        s = min_samples
        profiles = self.sort_profile()
        for m, field in profiles.items():
            time = field["Time"]
            daily = [t.split()[0] for t in time]   # Get days only
            days = list(set(daily))       # Different type of days
            days.sort()
            # List of tuple (index, day)
            inds = [(daily.index(d), d) for d in days
                    if (n_samples[m][d] >= s)]
            daily_i[m] = inds

        return daily_i

    def filter_daily_samples(self,
                             min_samples: int = 96) -> dict:
        """Filter uncompleted daily data.

        It returns profiles with at least ``min_samples``
        measures per day.

        """
        daily_i = self.get_iday(min_samples)
        n_samples = self.count_daily_samples()
        self_profiles = self._profiles     # Original profiles
        other_profiles = {}                     # Filtered profiles

        for m, inds in daily_i.items():
            times = []
            V_data = []
            P_data = []
            Q_data = []
            S_data = []
            for (i, day) in inds:
                j = i + n_samples[m][day]
                t = self_profiles[m]["Time"][i:j]
                v = self_profiles[m]["Voltage"][i:j]
                p = self_profiles[m]["Active"][i:j]
                q = self_profiles[m]["Reactive"][i:j]
                s = self_profiles[m]["Apparent"][i:j]
                times += t
                V_data += v
                P_data += p
                Q_data += q
                S_data += s

            other_profiles[m] = {
                "Time": times,
                "Voltage": V_data,
                "Active": P_data,
                "Reactive": Q_data,
                "Apparent": S_data
            }

        return other_profiles

    def daily_Vmax(self, min_samples: int = 85) -> dict[list]:
        """Maximum daily voltage per costumer.

        It gets the max voltage magnitud per load
        daily and the max active power registered at that instant.

        Note: In spite of the interval of time per sample
        is 15min the devices do not measure throughout the day
        so that only those with more than 90% of daily samples
        are considered hence ``samples=85``.

        """

        daily_i = self.get_iday(min_samples)
        n_samples = self.count_daily_samples()
        daily_max = {}
        for m, inds in daily_i.items():
            v = self._profiles[m]["Voltage"]
            p = self._profiles[m]["Active"]

            Vmax_t = []
            Pmax_t = []
            for (i, day_t) in inds:
                j = i + n_samples[m][day_t]
                if j == i:
                    continue
                else:
                    vt = v[i:j]
                    pt = p[i:j]
                Vmax = max(vt)
                i_max = vt.index(Vmax)
                Pmax = pt[i_max]

                Vmax_t.append(Vmax)
                Pmax_t.append(Pmax)
            daily_max[m] = {
                "Vmax": Vmax_t,
                "Pmax": Pmax_t
            }

        self._daily_Vmax = daily_max
        return daily_max

    def get_loadshapes(
            self, path: str = "./Perfiles/Puntos_de_medicion.xlsx",
            daily_samples: int = 96
    ) -> dict:
        """Generate loadshapes.

        It creats JSON file for each costumer whose demand
        satifies ``daily_samples`` per day while ``path``
        is the direcory of the file that contents the information
        about the transformer where the meter (load) is connected at.

        Note: It sets the attribute :py:attr:`Profile._loadshapes`.

        """
        samples = self.count_daily_samples()
        n = daily_samples
        loads = self.call_loads(path)
        replace_meterID = self.rename_meter()
        meters = [replace_meterID[M] for M in loads["Medidor"]]
        transformers = loads["Transformador"]
        tx_labels = list(set(transformers))
        profiles = self.sort_profile()
        ckt = {}          # Transformers ID as keys
        costumers = {}    # Meters as keys

        for m, data in profiles.items():
            daily = [time.split(" ")[0] for time in data["Time"]]
            days = list(set(daily))
            days.sort()
            # List of tuple (index, day)
            inds = [(daily.index(d), d) for d in days
                    if (samples[m][d] == n)]
            daily_curve = {}     # Days as keys
            for i, day in inds:
                j = i + samples[m][day]
                V = data["Voltage"][i:j]
                P = data["Active"][i:j]
                Q = data["Reactive"][i:j]
                S = data["Apparent"][i:j]
                curves = {
                    "Voltage": V,
                    "Active": P,
                    "Reactive": Q,
                    "Apparent": S
                }
                daily_curve[day] = curves

            costumers[m] = daily_curve

        for T in tx_labels:
            filter_meters = {}
            for c in costumers.keys():
                try:
                    m = replace_meterID[c]
                except KeyError as e:
                    print(f"NoCostumer:{e}")
                    continue
                tx_id = transformers[meters.index(m)]
                if T == tx_id:
                    filter_meters[m] = costumers[c]
            ckt[T] = filter_meters

        self._loadshapes = ckt
        return ckt

    def fill_vacuums(self) -> dict:
        """Deal with missing data.

        For those days missing in the data a demand curve
        of the closest day to it is allocated.

        """
        days = []
        ckt = self.get_loadshapes()
        # All possibles days
        for m, data in self._profiles.items():
            daily = [d.split(" ")[0] for d in data["Time"]]
            days += daily
        days = list(set(days))
        days.sort()

        for tx, m in ckt.items():
            for c, data in m.items():
                for day in days:
                    if day not in data:
                        d = 1
                        while (day not in data) and (d < 31):
                            delta = str(d)
                            t = pd.to_datetime(day, format='%d-%m-%Y')
                            dt = pd.Timedelta(f'{delta}d')
                            ti = t + dt
                            tj = t - dt
                            day_i = ti.strftime("%d-%m-%Y")
                            day_j = tj.strftime("%d-%m-%Y")
                            if day_i in data:
                                v = data[day_i]["Voltage"]
                                p = data[day_i]["Active"]
                                q = data[day_i]["Reactive"]
                                s = data[day_i]["Apparent"]
                                ckt[tx][c][day] = {
                                    "Voltage": v,
                                    "Active": p,
                                    "Reactive": q,
                                    "Apparent": s
                                }
                            elif day_j in data:
                                v = data[day_j]["Voltage"]
                                p = data[day_j]["Active"]
                                q = data[day_j]["Reactive"]
                                s = data[day_j]["Apparent"]
                                ckt[tx][c][day] = {
                                    "Voltage": v,
                                    "Active": p,
                                    "Reactive": q,
                                    "Apparent": s
                                }
                            else:
                                d += 1

                    else:
                        continue

        # Update attribute
        self._loadshapes = ckt
        return ckt

    def full_daily_Vmax(self) -> dict:
        ckt = self.fill_vacuums()
        daily_max = {}

        for tx, meters in ckt.items():
            for m, data in meters.items():
                Vmax_t = []
                Pmax_t = []
                for day in data.keys():
                    v = ckt[tx][m][day]["Voltage"]
                    Vmax = max(v)
                    Vmax_i = v.index(Vmax)
                    Pmax = ckt[tx][m][day]["Active"][Vmax_i]
                    Vmax_t.append(Vmax)
                    Pmax_t.append(Pmax)
                daily_max[m] = {
                    "Vmax": Vmax_t,
                    "Pmax": Pmax_t
                }
        # Update attribute
        self._daily_Vmax = daily_max
        return daily_max

    def set_daily_Edemand(self, freq: int = 4) -> dict:
        """Energy daily demand in [kWh].

        It considers samples are taken each 15min.
        It sets the field (column) ``Demand`` per day in [kWh]
        for all costumers in the data.

        Note: Multiply by (1/6) in case
        sample period is each 10min hence set freq parameter
        to six.

        """
        ckt = self.fill_vacuums()
        for tx, meters in ckt.items():
            for m, data in meters.items():
                for day in data.keys():
                    p = ckt[tx][m][day]["Active"]
                    kWh_day = 0
                    for kW in p:
                        # Interval each 15min
                        kWh_day += (kW/freq)
                    ckt[tx][m][day]["Demand"] = kWh_day

        return ckt

    def set_monthly_Edemand(self, n: int = 30) -> dict:
        """Energy monthly demand in [kWh].

        It sets firt the daily energy deman in [kWh]
        and then sets the field (column) ``KWHMONTH``
        for all costumers in the data. While ``n`` is the
        number of days in month to be considered.

        It sets the attribute :py:attr:`Profile._kWhmonth`.

        """
        ckt = self.set_daily_Edemand()

        demand_month = {}
        for tx, meters in ckt.items():
            for m, data in meters.items():
                demand_days = []
                t = []
                for d, day in enumerate(data.keys(), start=1):
                    if d <= n:
                        daily_demand = ckt[tx][m][day]["Demand"]
                        demand_days.append(daily_demand)
                        t.append(day)
                    else:
                        continue

                # Set day of max demand as the representative
                day_max = t[demand_days.index(max(demand_days))]
                v = ckt[tx][m][day_max]["Voltage"]
                p = ckt[tx][m][day_max]["Active"]
                q = ckt[tx][m][day_max]["Reactive"]
                month_kWh = round(sum(demand_days), 2)
                demand_month[m] = {
                    "KWHMONTH": month_kWh,
                    "Voltage": v,
                    "Active": p,
                    "Reactive": q
                }

        self._kWhmonth = demand_month
        return demand_month

    def set_curve_name_file(self):
        """Curve name attribute.

        It sets the attribute "CurveName" whose value
        would be the name of the file that contains
        the representative daily power damand, in the
        format bellow:

            "curve<kWhmonth>_<dd><class>"

        - kWhmonth: Energy monthly demand in kWh
        - dd: Two decimals
        - class: Sector (Commercial, Residential, Industry)

        """
        typical_curves = self.set_monthly_Edemand()
        for m, data in typical_curves.items():
            kWhmonth = round(data["KWHMONTH"], 2)
            name = str(kWhmonth)
            # To string ensuring the two decimals
            if ("." in name):
                if len(name.split(".")[1]) == 1:
                    name = f"{name}0"
            else:
                name = f"{name}.00"
            # File name attribute
            name = name.replace(".", "_")
            file_name = f"curve{name}{m[-1]}"
            file_name = file_name.lower()
            data["CurveName"] = file_name
            typical_curves[m] = data

        return typical_curves

    def write_curves(self,
                     dir: str = "./profiles/",
                     filef: str = "csv") -> None:
        """Representative daily curve.

        It creats ``filef`` (file format) files with
        the daily representative power
        within it.

        Note: Within the directory ``dir`` must be the
        empty folders: "residential", "industrial", "commercial".

        """
        curves = self.set_curve_name_file()
        for m, attrs in curves.items():
            cname = attrs["CurveName"]
            if "R" in m:
                path = f"{dir}residential/{cname}.{filef}"
            elif "C" in m:
                path = f"{dir}commercial/{cname}.{filef}"
            elif "I" in m:
                path = f"{dir}industrial/{cname}.{filef}"
            else:
                print(f"NoClass: {m}")
                continue
            # Write daily representative active power
            P = attrs["Active"]
            Q = attrs["Reactive"]
            with open(path, "w") as file:
                for p, q in zip(P, Q):
                    file.write(f"{p},{q}\n")


class Layers():

    def __init__(self):
        self._layers = []

    def get_shapes(self, path: str = "./Red/*.shp") -> None:
        """GIS data.

        It calls all shape files whose would be seen
        as layers and creates such attributes:

            1. Transformers
            2. Poles
            3. Meters
            4. Service Lines
            5. Secondary Low Voltage Lines.

        It also sets the same Coordinate Reference System (crs)
        at EPSG:5367 suitable for Costa Rica throughout
        layers.
        Note: The attribute :py:attr:`Layers._layers` contents
        :py:class:`geopandas.geodataframe.GeoDataFrame` object only.

        """
        files = glob.glob(path)
        for n, layer in enumerate(files):
            if n == 0:
                gdf01 = gpd.read_file(layer)
                self.transformers_gdf = gdf01.to_crs(epsg=5367)
                self._layers.append(self.transformers_gdf)
            elif n == 1:
                gdf02 = gpd.read_file(layer)
                self.poles_gdf = gdf02.to_crs(epsg=5367)
                self._layers.append(self.poles_gdf)
            elif n == 2:
                gdf03 = gpd.read_file(layer)
                self.meters_gdf = gdf03.to_crs(epsg=5367)
                self._layers.append(self.meters_gdf)
            elif n == 3:
                gdf04 = gpd.read_file(layer)
                self.services_lines_gdf = gdf04.to_crs(epsg=5367)
                self._layers.append(self.services_lines_gdf)
            elif n == 4:
                gdf05 = gpd.read_file(layer)
                self.LV_lines_gdf = gdf05.to_crs(epsg=5367)
                self._layers.append(self.LV_lines_gdf)

    def get_ckt_layers(
            self,
            profile: Profile,
            path: str = "../ModelBased/MAR/GIS/*.shp"
    ) -> None:
        shpfiles = glob.glob(path)
        # Remove duplicated layers
        shpfiles = [shp for shp in shpfiles if
                    ("Mario" not in shp)]
        shpfiles.sort()
        for (n, layer) in enumerate(shpfiles):
            # Aco_aerea
            if n == 0:
                gdf00 = gpd.read_file(layer)
                self.services_lines_gdf = gdf00.to_crs(epsg=5367)
                self._layers.append(self.services_lines_gdf)
            # Bus_BT_Layer
            elif n == 1:
                gdf01 = gpd.read_file(layer)
                self.LV_buses_gdf = gdf01.to_crs(epsg=5367)
                self._layers.append(self.LV_buses_gdf)
            # Bus_MT_Layer
            elif n == 2:
                gdf02 = gpd.read_file(layer)
                self.MV_buses_gdf = gdf02.to_crs(epsg=5367)
                self._layers.append(self.MV_buses_gdf)
            # Cargas
            elif n == 3:
                gdf03 = gpd.read_file(layer)
                self.LVloads_gdf_ckt = gdf03.to_crs(epsg=5367)
                self._layers.append(self.LVloads_gdf_ckt)
            # MT_aerea
            elif n == 4:
                gdf05 = gpd.read_file(layer)
                self.overH_MVlines_gdf = gdf05.to_crs(epsg=5367)
                self._layers.append(self.overH_MVlines_gdf)
            # MT_subterranea
            elif n == 5:
                gdf06 = gpd.read_file(layer)
                self.underG_MVlines_gdf = gdf06.to_crs(epsg=5367)
                self._layers.append(self.underG_MVlines_gdf)
            # Sec_aerea
            elif n == 6:
                gdf07 = gpd.read_file(layer)
                self.overH_LVlines_gdf = gdf07.to_crs(epsg=5367)
                self._layers.append(self.overH_LVlines_gdf)
            # Sec_subterranea
            elif n == 7:
                gdf08 = gpd.read_file(layer)
                self.underG_LVlines_gdf = gdf08.to_crs(epsg=5367)
                self._layers.append(self.underG_LVlines_gdf)
            # Subestacion
            elif n == 8:
                gdf09 = gpd.read_file(layer)
                self.subestation_gdf = gdf09.to_crs(epsg=5367)
                self._layers.append(self.subestation_gdf)
            # Transformadores
            elif n == 9:
                gdf10 = gpd.read_file(layer)
                self.transformers_gdf = gdf10.to_crs(epsg=5367)
                self._layers.append(self.transformers_gdf)

        # Profile loads shapes
        LVloads_gdf = profile.set_loads_attributes()
        self.LVloads_gdf = LVloads_gdf.to_crs(epsg=5367)
        self._layers.append(self.LVloads_gdf)

    def explore_map_ckt(self, profile: Profile) -> folium.Map:
        self.get_ckt_layers(profile)

        map = self.transformers_gdf.explore(
            tooltip=["PRIMCONN", "SECCONN"],
            popup=True,
            color="purple",
            name="transformers",
        )

        self.services_lines_gdf.explore(
            m=map,
            tooltip="TYPE",
            popup=True,
            color="aqua",
            name="serviceLines"
        )

        self.LV_buses_gdf.explore(
            m=map,
            tooltip="BUS",
            popup=True,
            color="red",
            name="LVbus"
        )

        self.MV_buses_gdf.explore(
            m=map,
            tooltip="BUS",
            popup=True,
            color="orange",
            name="MVbus"
        )

        self.LVloads_gdf_ckt.explore(
            m=map,
            tooltip="KWHMONTH",
            popup=True,
            color="magenta",
            name="LVloads_ckt"
        )

        self.LVloads_gdf.explore(
            m=map,
            tooltip="ICEobjID",
            popup=True,
            color="gold",
            name="LVloads"
        )

        self.overH_MVlines_gdf.explore(
            m=map,
            tooltip="NOMVOLT",
            popup=True,
            color="coral",
            name="overH_MVlines"
        )

        self.underG_MVlines_gdf.explore(
            m=map,
            tooltip="NOMVOLT",
            popup=True,
            color = "chocolate",
            name="underG_MVlines"
        )

        self.overH_LVlines_gdf.explore(
            m=map,
            tooltip="NOMVOLT",
            popup=True,
            color="maroon",
            name="overH_LVlines"
        )

        self.underG_LVlines_gdf.explore(
            m=map,
            tooltip="NOMVOLT",
            popup=True,
            color="orchid",
            name="underG_LVlines"
        )

        self.subestation_gdf.explore(
            m=map,
            tooltip="HIGHVOLT",
            popup=True,
            color="salmon",
            name="subestation"
        )

        folium.TileLayer("CartoDB positron", show=False).add_to(map)
        folium.LayerControl().add_to(map)

        self._map = map
        return map


class Graph():
    pass

    def split_groups():
        """Split apart Transformers.

        It creats sub-networks for each transformer.

        """
        pass


def set_service(n_phases: int) -> int:
    if n_phases == 2:
        return 12
    elif n_phases == 3:
        return 123
    else:
        print(f"NoService: {n_phases}")
        return None


def set_NOMVOLT(V: float) -> int:
    if V == 240.0:
        return 30
    else:
        print(f"NoVCode: {V}")
        return None


def explore_map(path: str = "./Red/*.shp") -> folium.Map:
    """Explore COOPELESCA network.

    Layers order
    0. AcometidasUCR
    1. MedidoresUCR
    2. PostesUCR
    3. SecundarioUCR
    4. TransformadoresUCR

    """
    shpfiles = glob.glob(path)
    shpfiles.sort()

    network = []
    for s in shpfiles:
        network.append(gpd.read_file(s))

    for n, layer in enumerate(network):
        # CR
        layer = layer.to_crs(epsg=5367)
        # Add layers to same folium
        if n == 0:
            map = layer.explore(
                name="serviceLines",
                color="salmon"
            )
        elif n == 1:
            layer.explore(
                m=map,
                name="meters",
                color="teal"
            )
        elif n == 2:
            layer.explore(
                m=map,
                name="LVbuses",
                color="violet"
            )
        elif n == 3:
            layer.explore(
                m=map,
                name="secondaryLines",
                color="tan"
            )
        elif n == 4:
            layer.explore(
                m=map,
                name="transformers",
                color="lime"
            )
    folium.TileLayer("CartoDB positron", show=False).add_to(map)
    folium.LayerControl().add_to(map)
    return map


if __name__ == "__main__":
    pass
