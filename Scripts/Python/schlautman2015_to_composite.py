#!/usr/bin/env python3

import argparse
import numpy as np
from datar import f
from datar.all import arrange, c, colnames, ungroup, rename, bind_rows
from datar.base import is_in
from datar.core.tibble import Tibble
from datar.dplyr import mutate, filter, inner_join, select, group_by, group_modify
from datar.tidyr.nest import nest,unnest
from datar.tibble import tibble
import pandas as pd
from scipy.interpolate import interp1d
import sys


class ssr2Composite():
    def _init_(self):
        return

    def setSSRMapPath(self, path):
        self.ssrmap = pd.read_csv(path, skiprows=1, index_col=0)
        self.ssrmap.columns = ["lg", "marker", "locus", "type", "position"]
        self.ssrmap = tibble(self.ssrmap) >> \
                        mutate(lg = lambda s: re.match(r"Cranberry-CNJ02-1-2007\.LG(\d+)", s).match(1),
                               position = lambda s: mm = re.match(r"([\d.]+)(\s*)(cM)*", s).match(1)) >> \
                        select("marker", "lg", "position")
        return

    def setCompositeMapPath(self, path):
        self.compositemap = tibble(pd.read_csv(path)) >> \
                                select("marker", "consensus", "LG") >> \
                                rename(position="consensus",
                                       lg="LG") >> \
                                select("marker", "lg", "position")
        return

    def genMap(self):
        self.combinedmap = self.ssrmap >> \
                            inner_join(self.compositemap, by="marker", suffix=("_ssr", "_composite"))
        self.map = {}
        for c in np.unique(self.combinedmap["lg"]):
            self.map[c] = interp1d(self.combinedmap >> filter(lg == c))

    def convert(self, lg, ssrCM):
        return(self.map[lg](ssrCM))
