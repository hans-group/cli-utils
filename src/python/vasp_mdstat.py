#!/usr/bin/env python
import re
import subprocess

import matplotlib.pyplot as plt
import numpy as np

plt.rc("font", size=13)

re_number_base = r"[\d+]?\+?-?.\d+?E[\+-]?\d+"
re_PE = re.compile(r"E0=\s+({})".format(re_number_base))
re_KE = re.compile(r"EK=\s+({})".format(re_number_base))
re_T = re.compile(r"T=\s+(\d+)")


def capture_float(re_, string):
    return np.fromiter(map(float, re_.findall(string)), dtype=np.float64)


def read_results():
    out = subprocess.check_output(["grep", "T=", "OSZICAR"]).decode("utf-8")

    T = capture_float(re_T, out)
    PE = capture_float(re_T, out)
    KE = capture_float(re_T, out)

    return dict(T=T, PE=PE, KE=KE)


if __name__ == "__main__":
    results = read_results()

    plt.figure(figsize=(6, 3.5))
    plt.xlabel("Step")
    plt.ylabel("Temperature (K)")
    plt.plot(results["T"])
    plt.tight_layout()
    plt.savefig("T.png", dpi=300)

    plt.figure(figsize=(6, 3.5))
    plt.xlabel("Step")
    plt.ylabel(r"E$_{pot}$ (eV)")
    plt.plot(results["PE"])
    plt.tight_layout()
    plt.savefig("PE.png", dpi=300)

    plt.figure(figsize=(6, 3.5))
    plt.xlabel("Step")
    plt.ylabel(r"E$_{kin}$ (eV)")
    plt.plot(results["KE"])
    plt.tight_layout()
    plt.savefig("KE.png", dpi=300)

    plt.figure(figsize=(6, 3.5))
    plt.xlabel("Step")
    plt.ylabel(r"E$_{tot}$ (eV)")
    plt.plot(results["KE"] + results["PE"])
    plt.tight_layout()
    plt.savefig("Etot.png", dpi=300)

    # Save data
    data = np.c_[results["T"], results["PE"], results["KE"]]
    np.savetxt("data.txt", data)
