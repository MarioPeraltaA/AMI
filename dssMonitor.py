"""Monitors data.

Parse, organize and store monitor data from OpenDSS.
"""


class dssData:
    def __init__(self):
        self.current_path = None
        self.dataPQ = []
        self.dataVI = []


    def monitor_PQ(self, PQpath: str, L: int) -> list[tuple]:
        """Get monitor power (P, Q) data.

        [kW] and [kVAr] where L is number of lines meassured by
        the monitor. It sums up the contribution of each phase.
        """
        self.current_path = PQpath
        route = PQpath
        monitorPQ = []
        with open(route, 'r') as f:
            for line in f:
                if 'hour' in line:
                    continue
                if len(line) > 1:
                    P = 0
                    Q = 0
                    data = line.split(',')
                    # Get channels
                    for p in range(1, L+1):
                        P += float(data[p*2].strip(' '))
                        Q += float(data[p*2 + 1].strip(' '))
                    monitorPQ.append((P, Q))
        self.dataPQ.append(monitorPQ)
        return monitorPQ


    def monitor_VI(self, VIpath: str, T: int) -> list[tuple[tuple]]:
        """Get monitor (V, A) data.

        [V] and [A] where T is number of threads meassured by
        the monitor. It sets the attribute dataVI as:
        [((Vabc, deg_v), (Iabc, deg_i))]
        """
        self.current_path = VIpath
        monitorVI = []
        route = VIpath
        with open(route, 'r') as f:
            for line in f:
                if 'hour' in line:
                    continue
                if len(line) > 1:
                    V_abc = []
                    I_abc = []
                    data = line.split(',')
                    # Get channels
                    for t in range(1, T+1):
                        v = float(data[2*t].strip(' '))
                        deg_v = float(data[2*t + 1].strip(' '))
                        c = 2*t + T*2
                        i = float(data[c].strip(' '))
                        deg_i = float(data[c+1].strip(' '))
                        V = (v, deg_v)
                        I = (i, deg_i)
                        V_abc.append(V)
                        I_abc.append(I)
                    monitorVI.append((V_abc, I_abc))

        self.dataVI.append(monitorVI)
        return monitorVI


if  __name__ == '__main__':
    pass
