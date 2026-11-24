"""Statistiques à l'échelle de l'unité expérimentale d'un banc de mesure."""

from pathlib import Path
from uuid import uuid4

import numpy as np
import pandas as pd
from scipy.stats import t


def sample(values, *, independent=True):
    """Moyenne, écart-type expérimental et intervalle de Student bilatéral.

    L'appelant fournit une observation par unité. Des appels répétés dans une
    seule exécution peuvent renseigner la dispersion conditionnelle. Dans ce
    cas, independent=False supprime l'intervalle entre unités indépendantes.
    """
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or not x.size or not np.isfinite(x).all():
        raise ValueError("Fournir une série unidimensionnelle finie et non vide.")
    count = len(x)
    mean = float(x.mean())
    sd = float(x.std(ddof=1)) if count > 1 else np.nan
    half = (float(t.ppf(.975, count - 1)) * sd / np.sqrt(count)
            if independent and count > 1 else np.nan)
    relative_half = half / abs(mean) if mean else np.nan
    return dict(n=count, mean=mean, sd=sd, low95=mean-half, high95=mean+half,
                relative_half95=relative_half,
                minimum=float(x.min()), maximum=float(x.max()))


class Measurements:
    """Conserver les résumés tracés, avec leur périmètre et leur effectif."""

    def __init__(self):
        self._records = {}

    def add(self, figure, metric, series, case, values, unit, *, independent=True):
        record = sample(values, independent=independent)
        key = (figure, metric, series, case)
        self._records[key] = dict(figure=figure, metric=metric, series=series, case=case,
                                  sampling_unit=unit, independent=independent, **record)
        return record

    def export(self, destination):
        data = pd.DataFrame(self._records.values())
        target = Path(destination)
        temporary = target.with_name(f".{target.name}.{uuid4().hex}.tmp")
        try:
            with temporary.open("x", encoding="utf-8", newline="") as stream:
                data.to_csv(stream, index=False, float_format="%.10g")
            temporary.replace(target)
        finally:
            temporary.unlink(missing_ok=True)
