using PyCall
using JLD2

py"""
import numpy as np
from batman import _rsky

sky_coords = {}
def sky_coords():
    t = np.linspace(-100, 100, 1_000)

    t0, period, a, e, omega, incl = (
        x.flatten()
        for x in np.meshgrid(
            np.linspace(-5.0, 5.0, 2),
            np.exp(np.linspace(np.log(5.0), np.log(50.0), 3)),
            np.linspace(50.0, 100.0, 2),
            np.linspace(0.0, 0.9, 5),
            np.linspace(-np.pi, np.pi, 3),
            np.arccos(np.linspace(0, 1, 5)[:-1]),
        )
    )

    r_batman = np.empty((len(t), len(t0)))

    for i in range(len(t0)):
        r_batman[:, i] = _rsky._rsky(
            t, t0[i], period[i], a[i], incl[i], e[i], omega[i], 1, 1
        )

    m = r_batman < 100.0

    return {
        "m_sum" : m.sum().item(), # Save native Int format
        "r_batman" : r_batman,
        "m" : m,
        "t" : t,
        "t0" : t0,
        "period" : period,
        "a" : a,
        "e" : e,
        "omega" : omega,
        "incl" : incl,
    }
"""
save("test_data/KeplerianOrbit_sky_coords.jld2", py"sky_coords"())
