"""Hosting Capacity ModelBase.

This program interacts with the API-Extention of
OpenDSS for Python in order to estimate the hosting
capacity based on simulated data of any clusters.

Such clusters are identify per transformer groups as follow:

    - Transformer TRA24198: Group1 (33 loads)
    - Transformer TRA21089: Group2 (26 loads)
    - Transformer TRA1624: Group3 (33 loads)

Author: Mario R. Peralta A.
email: Mario.Peralta@ucr.ac.cr

"""
from dss import DSS as DSSobj
import pandas as pd
import numpy as np
import json   # Read real data


def random_data(loadshape_name: str,
                N: int = 96,
                Q: bool = False) -> np.ndarray:
    """Artificial daily curve.

    Normal distribution if is residential load class
    and Uniform for commercial or industrial.
    In case random data for reactive power wants to be
    generated turn ``Q`` into True.
    Note: It assumes 96 meassures per day with interval
    of 15min between samples.

    """
    if Q:
        return None
    else:
        # Monthly energy demand kWh to Power
        kWhmonth_name = loadshape_name.replace("curve", "")[:-1]
        kWhmonth_val = float(kWhmonth_name.replace("_", "."))
        # Daily power
        dailyP = (kWhmonth_val*4) / 30
        # Amount per sample with 96 per day
        Pi = dailyP / N
        if loadshape_name[-1] == "r":
            Pdata = np.random.normal(Pi, 0.01, N)
        else:
            low = dailyP / N
            high = low + 0.001
            Pdata = np.random.uniform(low, high, N)

        return Pdata


def random_loadshapes(dssCircuit,
                      Q: bool = False):
    """Representative daily curve.

    Update the Pmult attribute of LoadShapes object
    for costumers of the active circuit.
    Turn Q into True to also generate random data
    for reactive power.

    """
    loadshape_names = dssCircuit.LoadShapes.AllNames
    # Set LoadShape object to first element
    _ = dssCircuit.LoadShapes.First   # Set it to "default"
    for _ in loadshape_names:
        i = dssCircuit.LoadShapes.Next
        if i:
            name_attr = dssCircuit.LoadShapes.Name
            dssCircuit.LoadShapes.Pmult = random_data(
                name_attr, Q=False
            )
    if Q:
        # Set back to the first element
        _ = dssCircuit.LoadShapes.First
        return None
    # Set back to the first element
    _ = dssCircuit.LoadShapes.First
    return dssCircuit


def get_cluster(data: dict, meterID: str) -> str:
    """Retrive the cluster.

    It gets the transformer key in order
    to index the representative curve
    a given day.

    """
    for T, meters in data.items():
        if meterID in meters:
            return T
    else:
        print("NoCluster")
        return None


def get_labels(loadshape_name: str,
               edemand_path: str = "./DSS/profiles/edemand.json",
               daily_path: str = "./DSS/profiles/daily.json") -> tuple[str]:

    with open(edemand_path, "r") as json_file:
        edemand = json.load(json_file)

    with open(daily_path, "r") as json_file:
        daily = json.load(json_file)

    for m, data in edemand.items():
        if data["CurveName"] == loadshape_name:
            cluster_label = get_cluster(daily, m)
            if cluster_label:
                return (daily, cluster_label, m)
            else:
                print(f"IslandMeter: {m}")
                return None


def get_loadshapes(
        date: str,
        loadshape_name: str,
        edemand_path: str = "./DSS/profiles/edemand.json",
        daily_path: str = "./DSS/profiles/daily.json",
        QkVAr: bool = True
):
    """Retrieve real measurements of power.

    The json file ``edemand.json`` has the loadshape_name
    attribute linked to a load (costumer) while ``daily.json``
    has the representative demand curves per day.

    Note: ``date`` arguments comes in format: "%d-%m-%Y"
    and ``loadshape_name`` is the dssCircuit.LoadShape.Name
    property.

    """
    labels = get_labels(loadshape_name,
                        edemand_path,
                        daily_path)
    if labels:
        curves_data, T, meter = labels
    else:
        print(f"KeyError: {loadshape_name} not found")
        return None
    if QkVAr:
        P = curves_data[T][meter][date]["Active"]
        Q = curves_data[T][meter][date]["Reactive"]
        return P, Q
    else:
        P = curves_data[T][meter][date]["Active"]
        return P


def set_loadshapes(dssCircuit,
                   date: str,
                   QkVAr: bool = True):
    """Representative daily curve.

    Update both Pmult & Qmult attribute of LoadShapes object
    for costumers of the active circuit.
    Turn QkVAr into False to modify Pmult only.

    """
    loadshape_names = dssCircuit.LoadShapes.AllNames
    # Set LoadShape object to first element
    _ = dssCircuit.LoadShapes.First   # Set it to "default"
    for _ in loadshape_names:
        name_attr = dssCircuit.LoadShapes.Name
        # Both P & Q
        if "curve" in name_attr:
            daily_cuves = get_loadshapes(date, name_attr)
            if daily_cuves:
                P, Q = daily_cuves
                # P & Q
                if QkVAr:
                    dssCircuit.LoadShapes.Pmult = P
                    dssCircuit.LoadShapes.Qmult = Q
                    _ = dssCircuit.LoadShapes.Next
                # P only
                else:
                    dssCircuit.LoadShapes.Pmult = P
                    _ = dssCircuit.LoadShapes.Next
            else:
                print(f"NoLoadShape: {name_attr}")
                return None
        else:
            _ = dssCircuit.LoadShapes.Next

    # Set back to the first element
    _ = dssCircuit.LoadShapes.First

    return dssCircuit


def set_VPQ_df(load_names, dti) -> tuple[pd.DataFrame]:
    # Define dataframe for each magnitud
    dataV_loads = pd.DataFrame(index=dti,
                               columns=load_names,
                               dtype=np.float64)
    dataP_loads = pd.DataFrame(index=dti,
                               columns=load_names,
                               dtype=np.float64)
    dataQ_loads = pd.DataFrame(index=dti,
                               columns=load_names,
                               dtype=np.float64)

    return (dataV_loads, dataP_loads, dataQ_loads)


def run_Daily(
        dssCircuit
) -> tuple[pd.DataFrame]:
    """Solve circuit.

    Run circuit each instant throught one day.
    Run a daily, and store its data in a DataFrame.
    Note: Assume single phase kind of loads only.

    """
    samples = 4 * 24   # Daily
    # Time domain
    dti = pd.date_range(start="00:00:00",
                        periods=samples,
                        freq="15T")
    dti = list(dti.strftime("%H:%M:%S"))

    load_names = dssCircuit.Loads.AllNames
    VPQ_df = set_VPQ_df(load_names, dti)
    dataV_loads = VPQ_df[0]
    dataP_loads = VPQ_df[1]
    dataQ_loads = VPQ_df[2]
    for t in dti:
        dssCircuit.Solution.Solve()
        for L in load_names:
            _ = dssCircuit.SetActiveElement(f"load.{L}")
            # Voltage
            load_element = dssCircuit.ActiveCktElement
            Va, _, Vb, _ = load_element.VoltagesMagAng
            V = Va + Vb
            # Power
            P, Q = load_element.TotalPowers
            # Update data
            dataV_loads.loc[t, L] = V
            dataP_loads.loc[t, L] = P
            dataQ_loads.loc[t, L] = Q

    return (dataV_loads, dataP_loads, dataQ_loads)


def split_groups(
        LVloads_path: str = "./DSS/MAR_LoadsLV.dss"
) -> dict[dict[list]]:
    """Data structure.

    Data structure to store constumers
    daily maximum values.
    It gets the path of *.dss file where the loads
    commands are at.

    """
    loads_dss = LVloads_path
    clusters = {}     # Transformers
    with open(loads_dss, "r") as file:
        for line in file:
            if "Load." in line:
                loadName = line.split(" ")[1].split(".")[1]
                loadName = loadName.lower()
                # Retrieve transformer group
                txlabel = line.split(" ")[-2]
                txlabel = txlabel.replace("=", "").strip("!")
                # Allocate load to transformer
                if txlabel in clusters:
                    clusters[txlabel].append(loadName)
                else:
                    clusters[txlabel] = [loadName]

    # Split into dicts the clusters
    ckt = {}
    for group, loads in clusters.items():
        costumer = {}
        for load in loads:
            # Daily max data of each costumer
            data_meter = {
                "Vmax": [],
                "Pmax": [],
                "Qmax": []
            }
            costumer[load] = data_meter
        ckt[group] = costumer

    return ckt


def get_Vmax(dssCircuit,
             Ndays: int = 30) -> tuple[list]:
    """Virtual data.

    Retrieve maximun voltage throughout one day
    and the power at that very instant and return
    as (Vmax, Pt, Qt) where each element of the tuple
    is a list of length equal to the number of ``Ndays``
    specified.

    """
    # Initialize data structure
    ckt = split_groups()
    # Run circuit during a month
    for a in range(1, Ndays+1):
        # From integer to date label
        day = f"{a}-07-2023"
        day_dti = pd.to_datetime(day, format="%d-%m-%Y")
        day = day_dti.strftime("%d-%m-%Y")
        dssCircuit = set_loadshapes(dssCircuit, day)
        V_data, P_data, _ = run_Daily(dssCircuit)

        # Retrieve daily data
        for load in V_data.columns:
            v_max = max(V_data[load].values)
            idx = V_data[load].idxmax()
            p_max = P_data.loc[idx, load]
            for group, L in ckt.items():
                if load in L:
                    ckt[group][load]["Vmax"].append(v_max)
                    ckt[group][load]["Pmax"].append(p_max)

    return ckt


def HC_cluster(
        dssCircuit,
        group: int = 3,
        Ndays: int = 30
) -> tuple[list]:
    """Hosting capacity per cluster.

    It runs ``Ndays`` daily, and returns the
    the data of maximun voltage per day againt
    the active power in that very instant.
    In case the LV group of transformer does not exist
    it will return ``None``.

    """
    cluster = f"Group{group}"
    ckt = split_groups()
    if cluster not in ckt:
        print(f"MissingCluster: {cluster}")
        return None
    else:
        hc = get_Vmax(dssCircuit, Ndays)
        hc_group = hc[cluster]
        vmax_data = []
        p_data = []
        for load in hc_group.values():
            vmax_data += load["Vmax"]
            p_data += load["Pmax"]

        return (p_data, vmax_data)


def get_LoadBuses(dssCircuit) -> list[tuple]:
    """Returns as (Load.name, LoadBusname).

    """
    LoadBuses = []
    for obj in dssCircuit.AllElementNames:
        if "Load" in obj:
            _ = dssCircuit.SetActiveElement(obj)
            load_at_bus = dssCircuit.ActiveCktElement.BusNames[0]
            LoadBuses.append((obj, load_at_bus))
        else:
            continue
    return LoadBuses


def creat_Vmonitor(loadName: str, dssText, dssCircuit):
    """Load to be sampled.

    Iterative monitor. Unique object the measure V of
    specific Load object. It should get reset each time
    moves to another load.

    """
    dssText.Command = f"new monitor.Vload Element={loadName} mode=0 ppolar=yes"
    dssCircuit.Solution.Solve()  # Non-PVSystem
    Va = dssCircuit.Monitors.AsMatrix()[:, 2]
    Vb = dssCircuit.Monitors.AsMatrix()[:, 4]
    V = Va + Vb
    return max(V)


def reset_monitor(loadName: str, dssCircuit):
    M = len(dssCircuit.Monitors.AllNames)
    _ = dssCircuit.Monitors.First
    for _ in range(M):
        if dssCircuit.Monitors.Name == "vload":
            dssCircuit.Monitors.Reset()
            dssCircuit.Monitors.Element = loadName
            dssCircuit.Solution.Solve()
            Va = dssCircuit.Monitors.AsMatrix()[:, 2]
            Vb = dssCircuit.Monitors.AsMatrix()[:, 4]
            V = Va + Vb
            return max(V)
        else:
            dssCircuit.Monitors.Next


def update_Vmax(dssCircuit):
    M = len(dssCircuit.Monitors.AllNames)
    _ = dssCircuit.Monitors.First
    for _ in range(M):
        if dssCircuit.Monitors.Name == "vload":
            Va = dssCircuit.Monitors.AsMatrix()[:, 2]
            Vb = dssCircuit.Monitors.AsMatrix()[:, 4]
            V = Va + Vb
            return max(V)
        else:
            dssCircuit.Monitors.Next


def redirect_PVsystem(dssText, loadBus: str) -> None:
    """Set ``PVtest`` object at first load bus.

    """
    dssText.Command = "redirect ./PVcurves.dss"
    dssText.Command = f"edit PVSystem.PVtest bus1={loadBus}"


def i_penetration(
        Vmax: float, constrain: float,
        dssText, dssCircuit
) -> tuple[float]:
    """Increase iteratively penetration.

    Stepsize of kVArated of += 1kVA for each itetarion.
    Fixed limit of 100 kVA.
    Retuns data as (Vmax, kVA).

    """
    # Already over constrain
    if Vmax > constrain:
        return (Vmax, 0)

    kVAs = range(1, 101)
    # Update penetration
    for s in kVAs:
        Vlim = Vmax
        Slim = s
        dssText.Command = f"edit PVSystem.PVtest kva={s}"
        dssCircuit.Solution.Solve()
        Vmax = update_Vmax(dssCircuit)
        if Vmax > constrain:
            return (Vlim, Slim)
        else:
            continue
    if Vmax < constrain:
        print(f"IterationLimitReached: {s}")
        return None


def c_penetration(
        Vmax, constrain: float,
        dssText, dssCircuit
) -> tuple[float]:
    """Increase iteratively penetration.

    Stepsize of kVArated of += 1kVA for each itetarion.
    Fixed limit of 100 kVA.
    Retuns data as (Vmax, kVA).

    Note: Updades loadshape (representative daily demand
    each iteration). Too slow.

    """
    # Already over constrain
    if Vmax > constrain:
        return (Vmax, 0)

    kVAs = range(1, 101)
    # Update penetration
    for a, s in enumerate(kVAs):
        if a <= 30:
            a += 1
        if 30 < a <= 60:
            a = a - 30
        if 60 < a <= 90:
            a = a - 60
        if 90 < a < 120:
            a = a - 90
        # From integer to date label
        day = f"{a}-07-2023"
        day_dti = pd.to_datetime(day, format="%d-%m-%Y")
        day = day_dti.strftime("%d-%m-%Y")
        dssCircuit = set_loadshapes(dssCircuit, day)
        if not dssCircuit:
            print("NoData")
            return None
        Vlim = Vmax
        Slim = s
        dssText.Command = f"edit PVSystem.PVtest kva={s}"
        dssCircuit.Solution.Solve()
        Vmax = update_Vmax(dssCircuit)
        if Vmax > constrain:
            return (Vlim, Slim)
        else:
            continue


def run_HC(dssText, dssCircuit) -> list[tuple]:
    LoadBuses = get_LoadBuses(DSScircuit)
    constrain = 240 * (1 + 0.05)
    maxData = []
    for n, load_obj in enumerate(LoadBuses):
        loadName = load_obj[0]
        loadBus = load_obj[1]
        # Set monitor Vload with no PV
        if n == 0:
            # Solve with zero penetration
            Vmax = creat_Vmonitor(loadName, dssText, dssCircuit)
            # Redirect PVSystem
            redirect_PVsystem(dssText, loadBus)
        else:
            Vmax = reset_monitor(loadName, dssCircuit)
            # Update PVbus
            dssText.Command = f"edit PVSystem.PVtest bus1={loadBus}"

        VSlim = i_penetration(Vmax,
                              constrain,
                              dssText,
                              dssCircuit)
        if VSlim:
            Vlim, Slim = VSlim
            meterID = loadName.split(".")[1]
            maxData.append((meterID, Vlim, Slim))
        else:
            maxData.append((VSlim))
            print(f"NoEnoughHC: {loadName}")
    return maxData


def allocate_HC(dssText, dssCircuit):
    """Store HC value to its cluster.

    """
    ckt = split_groups()
    maxData = run_HC(dssText, dssCircuit)
    for HCvals in maxData:
        meter, Vlim, Slim = HCvals
        cluster = get_cluster(ckt, meter)
        ckt[cluster][meter]["Vmax"] = Vlim
        ckt[cluster][meter]["Pmax"] = Slim

    return ckt


def write_HCtxt(dssText, dssCircuit) -> None:
    ckt = allocate_HC(dssText, dssCircuit)
    Tx = {
        "Group1": "TRA24198",
        "Group2": "TRA21089",
        "Group3": "TRA1624"
    }

    with open("./HC.txt", "w") as file:
        for group in ckt.keys():
            Tid = Tx[group]
            file.write("+--------------------------------------------+\n")
            file.write(f"|    Alojamiento en transformador {Tid}   |\n")
            file.write("+--------------------------------------------+\n")
            file.write(" No. &  Cliente  &   Vmax   &      Pmax       \\\ \hline \n")
            file.write("+--------------------------------------------+\n")
            for n, client in enumerate(ckt[group].keys(), start=1):
                meter = client
                Vmax = ckt[group][client]["Vmax"]
                Pmax = ckt[group][client]["Pmax"]
                file.write(f"   {n}  & ${meter}$ &  {Vmax:0.4f} &      {Pmax}  \\\ \n")
            file.write("+--------------------------------------------+\n")


def avgHC(dssText, dssCircuit) -> dict:
    """Avarage HC for each cluster.

    """
    ckt = allocate_HC(dssText, dssCircuit)
    Tx = {
        "Group1": "TRA24198",
        "Group2": "TRA21089",
        "Group3": "TRA1624"
    }
    HCs = {}
    for group, meters in ckt.items():
        HCvals = []
        cluster = Tx[group]
        for val in meters.values():
            HCvals.append(val["Pmax"])
        HCavg = sum(HCvals) / len(HCvals)
        HCs[cluster] = round(HCavg, 4)

    return HCs


def redirect_ckt(dssObj, stepsize: int = 15) -> tuple:
    """Call and set ckt.

    Where ``stepsize`` is the sample period
    in minutes for a daily.

    """
    # Call system
    dssText = dssObj.Text
    dssCircuit = dssObj.ActiveCircuit

    # Write new commands
    dssText.Command = "clear"
    dssText.Command = "New Circuit.MAR"
    dssText.Command = "set defaultbasefrequency=60"
    dssText.Command = "Edit Vsource.Source BasekV=138.0 pu=1.00"
    dssText.Command = "~ angle=0 frequency=60 phases=3"
    dssText.Command = 'Redirect DSS/MAR_OutputQGIS2OpenDSS.dss'
    dssText.Command = 'Set mode=daily'
    dssText.Command = 'Set number=1'
    dssText.Command = f'Set stepsize={stepsize}m'
    dssText.Command = 'Set time=(0,0)'
    # DSScircuit = set_loadshapes(DSScircuit, "06-07-2023")
    return dssText, dssCircuit


if __name__ == "__main__":
    DSStext, DSScircuit = redirect_ckt(DSSobj)
    # write_HCtxt(DSStext, DSScircuit)
    HCavg = avgHC(DSStext, DSScircuit)
    print(HCavg)
