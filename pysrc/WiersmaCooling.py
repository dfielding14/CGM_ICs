"""Float-native cooling tables used by the steady-state solver."""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np

import numpy_compat  # noqa: F401


BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "cooling" / "Wiersma09_CoolingTables"


def _broadcast_shape(*values):
    return np.broadcast(*[np.asarray(value, dtype=float) for value in values]).shape


def _safe_log(values):
    return np.log(np.clip(np.asarray(values, dtype=float), np.finfo(float).tiny, None))


def _closest_index(values: np.ndarray, target: float) -> int:
    return int(np.argmin(np.abs(np.asarray(values, dtype=float) - float(target))))


class Constant_Cooling:
    """Density-independent constant cooling function."""

    def __init__(self, lambda_cgs_value: float):
        self.lambda_cgs_value = float(lambda_cgs_value)

    def lambda_cgs(self, temperature_K, nH_cm3):
        shape = _broadcast_shape(temperature_K, nH_cm3)
        return np.full(shape, self.lambda_cgs_value, dtype=float)

    def dln_lambda_dln_T(self, temperature_K, nH_cm3):
        return np.zeros(_broadcast_shape(temperature_K, nH_cm3), dtype=float)

    def dln_lambda_dln_rho(self, temperature_K, nH_cm3):
        return np.zeros(_broadcast_shape(temperature_K, nH_cm3), dtype=float)


class Wiersma_Cooling:
    """Signed net cooling, bilinear in log T and log nH; clamp at table edges.

    Use Y=0.248, the primordial H/He table also shipped in GOTHAM. Metal rates
    scale linearly with gas-phase Z/Zsun, as in GOTHAM's optically thin limit.
    """

    def __init__(self, Z2Zsun: float, z: float):
        if not np.isfinite(Z2Zsun) or Z2Zsun < 0 or not np.isfinite(z):
            raise ValueError("metallicity must be nonnegative and redshift finite")
        files = []
        redshifts = []
        for table_path in sorted(DATA_DIR.glob("z_*.hdf5")):
            try:
                redshift = float(table_path.stem.removeprefix("z_"))
            except ValueError:
                continue
            redshifts.append(redshift)
            files.append(table_path)
        if not files:
            raise FileNotFoundError(f"no Wiersma cooling tables found under {DATA_DIR}")
        if not min(redshifts) <= z <= max(redshifts):
            raise ValueError("redshift is outside the Wiersma tables")
        index = _closest_index(redshifts, z)
        self.redshift = redshifts[index]
        with h5py.File(files[index], "r") as handle:
            helium_index = _closest_index(handle["Metal_free/Helium_mass_fraction_bins"][...], 0.248)
            self._rates = (np.asarray(handle["Metal_free/Net_Cooling"][helium_index, ...], dtype=float)
                           + Z2Zsun * np.asarray(handle["Total_Metals/Net_cooling"][...], dtype=float))
            self._logT = np.log(np.asarray(handle["Metal_free/Temperature_bins"][...], dtype=float))
            self._logn = np.log(np.asarray(handle["Metal_free/Hydrogen_density_bins"][...], dtype=float))

    @classmethod
    def from_config(cls, config):
        return cls(Z2Zsun=config.metallicity_solar, z=config.redshift)

    def _evaluate(self, temperature_K, nH_cm3, derivative=None):
        temperature, density = np.broadcast_arrays(np.asarray(temperature_K, dtype=float),
                                                  np.asarray(nH_cm3, dtype=float))
        if np.any(~np.isfinite(temperature) | ~np.isfinite(density) | (temperature <= 0) | (density <= 0)):
            raise ValueError("temperature and hydrogen density must be finite and positive")
        logT, logn = np.log(temperature), np.log(density)
        x, y = np.clip(logT, self._logT[0], self._logT[-1]), np.clip(logn, self._logn[0], self._logn[-1])
        i = np.clip(np.searchsorted(self._logT, x, side="right") - 1, 0, len(self._logT) - 2)
        j = np.clip(np.searchsorted(self._logn, y, side="right") - 1, 0, len(self._logn) - 2)
        dx, dy = self._logT[i + 1] - self._logT[i], self._logn[j + 1] - self._logn[j]
        t, n = (x - self._logT[i]) / dx, (y - self._logn[j]) / dy
        f00, f10 = self._rates[i, j], self._rates[i + 1, j]
        f01, f11 = self._rates[i, j + 1], self._rates[i + 1, j + 1]
        rate = (1 - t) * ((1 - n) * f00 + n * f01) + t * ((1 - n) * f10 + n * f11)
        if derivative is None:
            return rate
        if np.any(rate <= 0):
            raise ValueError("logarithmic cooling slopes require positive net cooling")
        if derivative == 0:
            slope = ((1 - n) * (f10 - f00) + n * (f11 - f01)) / (dx * rate)
            return np.where((logT < self._logT[0]) | (logT > self._logT[-1]), 0.0, slope)
        slope = ((1 - t) * (f01 - f00) + t * (f11 - f10)) / (dy * rate)
        return np.where((logn < self._logn[0]) | (logn > self._logn[-1]), 0.0, slope)

    def lambda_cgs(self, temperature_K, nH_cm3):
        return self._evaluate(temperature_K, nH_cm3)

    def dln_lambda_dln_T(self, temperature_K, nH_cm3):
        return self._evaluate(temperature_K, nH_cm3, derivative=0)

    def dln_lambda_dln_rho(self, temperature_K, nH_cm3):
        return self._evaluate(temperature_K, nH_cm3, derivative=1)


class Kartick_Cooling:
    """One-dimensional CIE cooling table used for comparison runs."""

    table_path = BASE_DIR / "cooling" / "Kartick_CIE_cooling.table"

    def __init__(self):
        table = np.genfromtxt(self.table_path, dtype=float)
        self.temperature_bins_K = table[:, 0]
        # Table stores Lambda * n_e. The original code converts this to the
        # n_H^2 Lambda convention used throughout the solver.
        self.lambda_bins_cgs = table[:, 1] * table[:, 2]

        log_temperature = _safe_log(self.temperature_bins_K)
        dlogT = np.diff(log_temperature)
        log_lambda = np.log(np.clip(self.lambda_bins_cgs, np.finfo(float).tiny, None))
        self._gradient_temperature = 0.5 * (log_temperature[1:] + log_temperature[:-1])
        self._gradient_values = np.diff(log_lambda) / dlogT

    def lambda_cgs(self, temperature_K, nH_cm3=None):
        temperature = np.asarray(temperature_K, dtype=float)
        values = np.interp(
            _safe_log(temperature),
            _safe_log(self.temperature_bins_K),
            np.log(np.clip(self.lambda_bins_cgs, np.finfo(float).tiny, None)),
        )
        return np.exp(values)

    def dln_lambda_dln_T(self, temperature_K, nH_cm3=None):
        return np.interp(_safe_log(temperature_K), self._gradient_temperature, self._gradient_values)

    def dln_lambda_dln_rho(self, temperature_K, nH_cm3=None):
        return np.zeros(np.shape(np.asarray(temperature_K, dtype=float)), dtype=float)


class DopitaSutherland_CIE:
    """Dopita & Sutherland (1996) CIE cooling table."""

    table_path = BASE_DIR / "cooling" / "DopitaSutherland_CIE.dat"

    def __init__(self, Z2Zsun: float):
        table = np.genfromtxt(self.table_path, dtype=float)
        self.temperature_bins_K = table[:, 0]
        if Z2Zsun == 1.0:
            lambda_bins = table[:, 1]
        elif Z2Zsun == 1.0 / 3.0:
            lambda_bins = table[:, 2]
        else:
            raise ValueError("DopitaSutherland_CIE only supports Z/Zsun = 1 or 1/3")

        n_i_to_n_H = 1.22 ** -1 / 0.7
        n_e_to_n_H = 1.17 ** -1 / 0.7
        self.lambda_bins_cgs = lambda_bins * n_i_to_n_H * n_e_to_n_H

        log_temperature = _safe_log(self.temperature_bins_K)
        dlogT = np.diff(log_temperature)
        log_lambda = np.log(np.clip(self.lambda_bins_cgs, np.finfo(float).tiny, None))
        self._gradient_temperature = 0.5 * (log_temperature[1:] + log_temperature[:-1])
        self._gradient_values = np.diff(log_lambda) / dlogT

    def lambda_cgs(self, temperature_K, nH_cm3=None):
        return np.exp(
            np.interp(
                _safe_log(temperature_K),
                _safe_log(self.temperature_bins_K),
                np.log(np.clip(self.lambda_bins_cgs, np.finfo(float).tiny, None)),
            )
        )

    def dln_lambda_dln_T(self, temperature_K, nH_cm3=None):
        return np.interp(_safe_log(temperature_K), self._gradient_temperature, self._gradient_values)

    def dln_lambda_dln_rho(self, temperature_K, nH_cm3=None):
        return np.zeros(np.shape(np.asarray(temperature_K, dtype=float)), dtype=float)
