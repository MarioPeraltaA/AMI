"""OpenDSS COM.

Run OpenDss from python

    Author: Mario Roberto Peralta A.
    email: Mario.Peralta@ucr.ac.cr

"""
from dss import DSS


dssText = DSS.Text

# # Compile test PVckt: 3loads
# dssText.Command = "compile ./dssfilesT/PVtest_ckt.dss"

# Compile 14loads PVckt
dssText.Command = "compile ./dssfiles/PV_ckt.dss"
