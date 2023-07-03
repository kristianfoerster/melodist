import unittest

import numpy as np
import pandas as pd

import melodist


class Test5MinutePrecipitation(unittest.TestCase):
    def test_cascade(self):
        data = pd.read_csv(
            'examples/testdata_precip5min.csv.gz',
            names=['time', 'precip'],
            parse_dates=True,
            index_col=0,
        )
        precip_daily = data.precip.resample('D').sum()

        cascopt = melodist.build_casc(data, hourly=False, level=9, percentile=90)
        precip_disagg = melodist.disagg_prec_cascade(
            precip_daily,
            cascopt[0],
            hourly=False,
            level=9,
        )

        assert np.allclose(
            precip_daily,
            precip_disagg.resample('D').sum(),
            atol=1e-3,
            equal_nan=True,
        )
