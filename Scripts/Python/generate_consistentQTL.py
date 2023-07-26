#!/usr/bin/env python3

#Merge common QTL across years

import operator
import pandas as pd
import numpy as np
import re
from datar import f
from datar.all import arrange, c, colnames, ungroup, rename, bind_rows
from datar.base import is_in, all_
from datar.base import any_
from datar.core.tibble import Tibble
from datar.dplyr import mutate, filter, inner_join, select, group_by, group_modify
from datar.tidyr.nest import nest,unnest
from datar.tibble import tibble
from pipda import register_verb, Context
import networkx as nx
from networkx.algorithms.clique import find_cliques
from networkx.algorithms.approximation import clique_removal

from datar.core.backends.pandas import DataFrame, Series
from datar.core.contexts import Context
from datar.core.utils import logger
from datar.core.tibble import Tibble, TibbleGrouped, reconstruct_tibble
from datar.core.broadcast import broadcast_to
from datar.core.operator import _binop

from scipy.interpolate import interp1d

columns_highlighted = ("method", "trait", "chr", "model", "position", "nearest_marker", "qtl_lod", "marker_variance", "model_variance", "interval_left", "interval_right", "population","study")
columns_highlighted_i = tuple(["id"] + list(columns_highlighted))
columns_effects_extended = ("effect_mean_AD","effect_mean_AC","effect_mean_BC","effect_mean_BD","effect_se_AD","effect_se_AC","effect_se_BC","effect_se_BD","AvB","CvD","Int")
columns_output = ("id", "trait", "model", "chr", "position", "mean_marker_variance", "mean_model_variance", "mean_combined_variance", "population", "method", "study", "interval_left", "interval_right", "count", "model_count", "study_count", "trait_count", "method_count", "population_count")
columns_output_extended = tuple(["mean_" + o for o in columns_effects_extended])

@register_verb(types=Tibble,context=Context.EVAL)
def group_trait(data,**kwargs):
    property_trait = kwargs['property_trait']
    property_name = kwargs['name']
    groups = kwargs['groups']
    traits = data[property_trait]
    traits_i = pd.Index(traits).factorize()
    traits_i_size = np.max(traits_i[0])+1
    group_categories = pd.Series(traits_i[0])
    for i,grouped_traits in enumerate(groups):
        group_categories[is_in(traits,grouped_traits)] = i+traits_i_size #Added np.where b/c indexes may not be contiguous or in range due to filtering on input dataset
    data =data.copy()
    data.index = group_categories.index
    data[property_name] = group_categories
    return(data)

@register_verb(types=Tibble,context=Context.EVAL)
def long2short(data,**kwargs):
    property_trait = kwargs['trait']
    property_map = kwargs['l2smap']
    data.index = range(0,data.shape[0])
    long = data[property_trait]
    short = np.ones(shape=(long.shape[0]),dtype="object")
    for i,l in enumerate(long):
        short[i] = property_map[l]
    data[property_trait] = pd.Series(short)
    return(data)

algo = 1
def merge_qtl(df):
    fstack = df >> arrange(f.interval_left)
    nstack = []
    goverlap = nx.Graph()
    for i in range(0,len(fstack)):
        goverlap.add_node(i,**(fstack.iloc[i,].to_dict()))
    for i in range(0,len(fstack)):
        qc = fstack.iloc[i,]
        for j in range(i-1,-1,-1):
            qp = fstack.iloc[j,]
            if (qc.interval_left < qp.interval_right):
                goverlap.add_edge(i,j)
    #Now that we've generated base graph, find all cliques, which represent the representative overlaps
    if algo == 1:
        _,maximal_cliques = clique_removal(goverlap)
    else:
        maximal_cliques  = find_cliques(goverlap)
    for clique in maximal_cliques:
        intervals_left = [goverlap.nodes[i]["interval_left"] for i in clique]
        intervals_right= [goverlap.nodes[i]["interval_right"] for i in clique]
        merger = dict((k,set([])) for k in goverlap.nodes[list(clique)[0]].keys())
        merger["marker_variance"] = [] #Marker variance needs to be a list and not a set, in case of repeated variance values
        merger["model_variance"] = [] #Model variance needs to be a list and not a set, in case of repeated variance values
        merger["combined_variance"] = []
        merger["interval_left"] = np.max(intervals_left)
        merger["interval_right"] = np.min(intervals_right)
        merger['position'] = (merger['interval_left'] + merger['interval_right'])/2
        if all_(is_in(columns_effects_extended,list(merger.keys()))):
            for k in columns_effects_extended:
                merger[k] = []
        for n in [goverlap.nodes[i] for i in clique]:
            n = dict([(k,set([v])) for k,v in n.items()]) #Convert to set for merger
            for k in n.keys():
                if( k == "marker_variance" ):
                    merger[k].append(list(n[k])[0])
                elif( k == "model_variance" ):
                    if not np.isnan(list(n[k])[0]):
                        merger[k].append(list(n[k])[0])
                        merger["combined_variance"].append((float(list(n["marker_variance"])[0]) * float(list(n[k])[0]))/100) #divide by 100 b/c each factor is a percent
                elif (k in columns_effects_extended):
                    if not np.isnan(list(n[k])[0]):
                        merger[k].append(list(n[k])[0])
                elif( k not in ['interval_left', 'interval_right', 'position'] ):
                    merger[k] = merger[k] | n[k] #Intersect sets
        #Convert all sets to string-concatenated lists
        for k in merger.keys():
            if( k not in ['interval_left', 'interval_right', 'position'] ):
                merger[k] = '+'.join([str(e) for e in merger[k]])
        merger['count'] = len(clique)
        merger['mean_marker_variance'] = (eval(merger['marker_variance']))/merger['count']
        if merger['model_variance'] != "":
            merger['mean_model_variance'] = (eval(merger['model_variance']))/merger['count']
            merger['mean_combined_variance'] = (eval(merger['combined_variance']))/merger['count']
        else:
            merger['mean_model_variance'] = np.nan
            merger['mean_combined_variance'] = np.nan
        merger['model_count'] = len(merger['model'].split('+'))
        merger['study_count'] = len(merger['study'].split('+'))
        merger['method_count'] = len(merger['method'].split('+'))
        merger['population_count'] = len(merger['population'].split('+'))
        merger['trait_count'] = len(merger['trait'].split('+'))
        if all_(is_in(columns_effects_extended,list(merger.keys()))):
            for k in columns_effects_extended:
                try:
                    merger['mean_'+k] = (eval(merger[k])/merger['count'])
                except SyntaxError as se:
                    print("Invalid syntax for key {}: \"{}\".  Merger: \"{}\"".format(k,merger[k],merger))
        nstack.append(merger)
    nstack = tibble(pd.DataFrame(nstack))
    return(nstack)

@register_verb(types=Tibble,context=Context.EVAL)
def filter_qtl(_data, *filters):
    _data = _data.copy()
    if _data.shape[0] != 0 and filters:
        condition = np.array(True)
        for cond in filters:
            condition = _binop(operator.and_, condition, cond)
        #for f in filters:
        #    _data = _data[f]
        grouper = None
        if isinstance(_data, TibbleGrouped):
            grouper = _data._datar["grouped"].grouper
        condition = broadcast_to(condition, _data.index, grouper)
        if isinstance(condition, np.bool_):
            condition = bool(condition)
        if condition is True:
            return _data.copy()
        if condition is False:
            return _data.take([])
        if isinstance(condition, Series):
            condition = condition.values
        _data = ungroup(_data, __ast_fallback="normal")[condition]
    if all_(is_in(columns_output_extended,list(_data.keys()))):
        co = list(columns_output) + list(columns_output_extended)
    else:
        co = list(columns_output)
    _data = _data[list(co)]
    _data.sort_values(by="mean_marker_variance", ascending=False, inplace=True, ignore_index=True)
    return(_data)

@register_verb(types=Tibble,context=Context.EVAL)
def fix_ssr_lg(_data, *kwargs):
    _data = _data.copy()
    def ssr_lg_replace(lg):
        return(re.match(r"Cranberry-CNJ02-1-2007\.LG(\d+)", lg)[1])
    vssr_lg_replace = np.vectorize(ssr_lg_replace)
    if _data.shape[0] != 0:
        _data.lg = vssr_lg_replace(_data.lg)
    return(_data)

@register_verb(types=Tibble,context=Context.EVAL)
def fix_ssr_position(_data, *kwargs):
    _data = _data.copy()
    def ssr_position_replace(position):
        return(re.match(r"([\d.]+)(\s*)(cM)*", position)[1])
    vssr_position_replace = np.vectorize(ssr_position_replace)
    if _data.shape[0] != 0:
        _data.position = vssr_position_replace(_data.position)
    return(_data)

@register_verb(types=Tibble,context=Context.EVAL)
def convert_ssr2composite(_data, **kwargs):
    ssr2c = kwargs['ssr2c_class']
    _data = _data.copy()
    def ssr2comp(lg, position_ssr):
        return(ssr2c.convert(lg, position_ssr))
    vssr2comp = np.vectorize(ssr2comp)
    if _data.shape[0] != 0:
        _data.position = vssr2comp(_data.chr, _data.position)
        _data.interval_left = vssr2comp(_data.chr, _data.interval_left)
        _data.interval_right = vssr2comp(_data.chr, _data.interval_right)
    return(_data)

class ssr2Composite():
    def _init_(self):
        return

    def setSSRMapPath(self, path):
        self.ssrmap = pd.read_csv(path, skiprows=1, index_col=0)
        self.ssrmap.columns = ["lg", "marker", "locus", "type", "position"]
        self.ssrmap = tibble(self.ssrmap) >> \
                        fix_ssr_lg() >> \
                        fix_ssr_position() >> \
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
        for c in np.unique(self.combinedmap.lg_ssr):
            curr = self.combinedmap >> filter(f.lg_ssr == c)
            self.map[int(c)] = interp1d(curr.position_ssr.astype("float"), curr.position_composite.astype("float"))
        return

    def convert(self, lg, ssrCM):
        return(self.map[lg](ssrCM))

ssr2c = ssr2Composite()
ssr2c.setSSRMapPath("../../Data/genetic_data/RawData/schlautman2015_ssr_map.csv")
ssr2c.setCompositeMapPath("../../Data/genetic_data/RawData/consensusMapAll2.csv")
ssr2c.genMap()

#Read in model config file that maps trait names to abbreviated labels
model_config = tibble(pd.read_csv("../../Data/publication/tables/model-traits.cfg.csv"))
model_traits = pd.Index(model_config.trait).factorize()
model_short = pd.Index(model_config.label_short).factorize()
model_trait_map = dict(zip(model_traits[1],model_short[1]))
interval_max = 15 #Eliminates wide QTL whose confidence interval is too large
mvariance_min = 10 #Eliminates QTL with smaller marker variances


#My mapping CNJ02
cnj02_groups=[["MFM","UMFM","UTBM","UBW","UBM"],["TY","SFY","PAC"],["UBL","UKLvW","UKEC","ULvW"],["UL","UDM"],["UNP","UNAF"]] #Group QTL by traits that have similar meaning/correlation
qtls_cnj02 = tibble(pd.read_csv("../../Workflows/9/traits/effects_collated.wide.csv")) >> \
                filter((f.chr2 == '') | np.isnan(f.chr2)) >> \
                filter((f.interval <= interval_max) & (f.marker_variance >= mvariance_min)) >> \
                mutate(population="cnj02") >> \
                mutate(study="self") >>  \
                select(columns_highlighted+columns_effects_extended) >> \
                long2short(trait='trait',l2smap=model_trait_map)

#nis = np.min(qtls_cnj02.index) #next index start
#nie = np.max(qtls_cnj02.index) + 1 #next index end
nis = 0
nie = qtls_cnj02.shape[0]
qtls_cnj02 = qtls_cnj02 >> \
                mutate(id=range(nis,nie))
qtls_cnj02.index = range(0,len(qtls_cnj02))
qtls_cnj02.to_csv("../../Data/publication/tables/cnj02_qtl_collated.raw.consensus.csv", index=False)

qtls_s1 = qtls_cnj02 >> \
            filter(f.method == "scanone") >> \
            group_trait(property_trait='trait', name='group_category', groups=[]) >> \
            group_by(f.group_category,f.chr) >> \
            group_modify(merge_qtl) >> \
            ungroup() >> \
            filter_qtl()

qtls_s1.to_csv("../../Data/publication/tables/cnj02_qtl_collated.s1.consensus.csv", index=False)

qtls_sw = qtls_cnj02 >> \
            filter(f.method == "stepwiseqtl") >> \
            group_trait(property_trait='trait', name='group_category', groups=[]) >> \
            group_by(f.group_category,f.chr) >> \
            group_modify(merge_qtl) >> \
            ungroup() >> \
            filter_qtl()

qtls_sw.to_csv("../../Data/publication/tables/cnj02_qtl_collated.sw.consensus.csv", index=False)

qtls_m = qtls_cnj02 >> \
            group_trait(property_trait='trait', name='group_category', groups=[]) >> \
            group_by(f.group_category,f.chr) >> \
            group_modify(merge_qtl) >> \
            ungroup() >> \
            filter_qtl()

qtls_m.to_csv("../../Data/publication/tables/cnj02_qtl_collated.consensus.csv", index=False)

qtls_nay = qtls_cnj02 >> \
            filter(f.model != 'all-years') >> \
            group_trait(property_trait='trait', name='group_category', groups=[]) >> \
            gropby(f.group_category,f.chr) >> \
            group_modify(merge_qtl) >> \
            ungroup() >> \
            filter_qtl()

qtls_nay.to_csv("../../Data/publication/tables/cnj02_qtl_collated.consensus.nay.csv", index=False)

qtls_g = qtls_cnj02 >> \
            group_trait(property_trait='trait', name='group_category', groups=cnj02_groups) >> \
            group_by(f.group_category,f.chr) >> \
            group_modify(merge_qtl) >> \
            ungroup() >> \
            filter_qtl()
#            filter_qtl(f.count >= 4, f.mean_marker_variance > 5)

qtls_g.to_csv("../../Data/publication/tables/cnj02_qtl_collated.grouped.consensus.csv", index=False)

#My mapping CNJ04
cnj04_groups=[["UNP","UNAF"],["USL","UDM"],["UNS","UBW","MFM","UBM""UMFM"],["TY","SFY","UNB","UTBM"],["UBL","UKLvW","UKEC","ULvW"]] #Group QTL by traits that have similar meaning/correlation


qtls_cnj04 = tibble(pd.read_csv("../../Workflows/10/traits/effects_collated.wide.csv")) >> \
                filter((f.chr2 == '') | np.isnan(f.chr2)) >> \
                filter((f.interval <= interval_max) & (f.marker_variance >= mvariance_min)) >> \
                mutate(population="cnj04") >> \
                mutate(study="self") >> \
                select(columns_highlighted+columns_effects_extended) >> \
                long2short(trait='trait',l2smap=model_trait_map)

nis = nie #Previous value
nie = nis + qtls_cnj04.shape[0]
qtls_cnj04 = qtls_cnj04 >> \
                mutate(id=range(nis,nie))
qtls_cnj04.index = range(0,len(qtls_cnj04))
qtls_cnj04.to_csv("../../Data/publication/tables/cnj04_qtl_collated.raw.consensus.csv", index=False)

qtls_s1 = qtls_cnj04 >> \
            filter(f.method == "scanone") >> \
            group_trait(property_trait='trait', name='group_category', groups=[]) >> \
            group_by(f.group_category,f.chr) >> \
            group_modify(merge_qtl) >> \
            ungroup() >> \
            filter_qtl()

qtls_s1.to_csv("../../Data/publication/tables/cnj04_qtl_collated.s1.consensus.csv", index=False)

qtls_sw = qtls_cnj04 >> \
            filter(f.method == "stepwiseqtl") >> \
            group_trait(property_trait='trait', name='group_category', groups=[]) >> \
            group_by(f.group_category,f.chr) >> \
            group_modify(merge_qtl) >> \
            ungroup() >> \
            filter_qtl()

qtls_sw.to_csv("../../Data/publication/tables/cnj04_qtl_collated.sw.consensus.csv", index=False)

qtls_m = qtls_cnj04 >> \
            group_trait(property_trait='trait', name='group_category', groups=[]) >> \
            group_by(f.group_category,f.chr) >> \
            group_modify(merge_qtl) >> \
            ungroup() >> \
            filter_qtl()
#            filter_qtl(f.count >= 2, f.mean_marker_variance > 5)

qtls_m.to_csv("../../Data/publication/tables/cnj04_qtl_collated.consensus.csv", index=False)

qtls_nay = qtls_cnj04 >> \
            filter(f.model != 'all-years') >> \
            group_trait(property_trait='trait', name='group_category', groups=[]) >> \
            group_by(f.group_category,f.chr) >> \
            group_modify(merge_qtl) >> \
            ungroup() >> \
            filter_qtl()
#            filter_qtl(f.count >= 2, f.mean_marker_variance > 5)

qtls_nay.to_csv("../../Data/publication/tables/cnj04_qtl_collated.consensus.nay.csv", index=False)

qtls_g   = qtls_cnj04 >> \
                group_trait(property_trait='trait', name='group_category', groups=cnj04_groups) >> \
                group_by(f.group_category,f.chr) >> \
                group_modify(merge_qtl) >> \
                ungroup() >> \
                filter_qtl()
#                filter_qtl(f.count >= 4, f.mean_marker_variance > 5)

qtls_g.to_csv("../../Data/publication/tables/cnj04_qtl_collated.grouped.consensus.csv", index=False)

#Merged/combined cnj02 and cnj04
cnj0x_groups=[["UNP","UNAF"],["MFM","UMFM","UBW","UBM"],["TY","SFY"],["UBL","UKLvW","UKEC","ULvW"]]
qtls_g   = qtls_cnj02 >> \
                bind_rows(qtls_cnj04) >> \
                group_trait(property_trait='trait', name='group_category', groups=cnj0x_groups) >> \
                group_by(f.group_category,f.chr) >> \
                group_modify(merge_qtl) >> \
                ungroup() >> \
                filter_qtl()
#                filter_qtl(f.count >= 4, f.mean_marker_variance > 5)

qtls_g.to_csv("../../Data/publication/tables/cnj0x_qtl_collated.grouped.consensus.csv", index=False)

#Luis Massive Phenotyping Dataset
mvariance_min = 8
massive_groups=[["BCOLOR","BCOLORVAR","TACY_DIFF","TACY_OCT","TACY_SEP"]] #Group QTL by traits that have similar meaning/correlation
qtls = tibble(pd.read_csv("../../Data/publication/tables/diazGarcia2018MassivePhenotyping.supplemental.qtls.csv"))
qtls_massive = qtls >> \
                    mutate(interval_left=f.position-(f.lod_15/2),
                           interval_right=f.position+(f.lod_15/2),
                           model_variance=np.nan,
                           qtl_lod=np.nan) >> \
                    filter((f.lod_15 <= interval_max) & (f.marker_variance >= mvariance_min)) >> \
                    mutate(study="diazGarcia2018_massive") >> \
                    select(columns_highlighted) >> \
                    mutate(count = 1)

nis = nie #Previous value
nie = nis + qtls_massive.shape[0]
qtls_massive = qtls_massive >> \
                mutate(id=range(nis,nie))
qtls_massive.index = range(0,len(qtls_massive))
qtls_massive.to_csv("../../Data/publication/tables/diazGarcia2018MassivePhenotyping.supplemental.raw.consensus.qtls.csv", index=False)

qtls_g = qtls_massive >> \
            group_trait(property_trait='trait', name='group_category', groups=massive_groups) >> \
            group_by(f.group_category,f.chr) >> \
            group_modify(merge_qtl) >> \
            ungroup() >> \
            filter_qtl()
#            filter_qtl(f.count >= 2, f.mean_marker_variance > 10.0)

qtls_g.to_csv("../../Data/publication/tables/diazGarcia2018MassivePhenotyping.supplemental.grouped.consensus.qtls.csv", index=False)


#Luis Image Phenotyping Dataset
mvariance_min = 6
image_groups=[["EC","LvW"]] #Group QTL by traits that have similar meaning/correlation
qtls = tibble(pd.read_csv("../../Data/publication/tables/diazGarcia2018ImagePhenotyping.supplemental.qtls.csv"))
qtls_image = qtls >> \
                mutate(interval_left=f.position-(f.lod_15/2),
                       interval_right=f.position+(f.lod_15/2),
                       model_variance=np.nan,
                       qtl_lod=np.nan,
                       population="gryg",
                       method="stepwiseqtl") >> \
                filter((f.lod_15  <= interval_max) & (f.marker_variance >= mvariance_min)) >> \
                mutate(study="diazGarcia2018_image") >> \
                select(columns_highlighted) >> \
                mutate(count = 1)

nis = nie #Previous value
nie = nis + qtls_image.shape[0]
qtls_image = qtls_image >> \
                mutate(id=range(nis,nie))
qtls_image.index = range(0,len(qtls_image))
qtls_image.to_csv("../../Data/publication/tables/diazGarcia2018ImagePhenotyping.supplemental.raw.consensus.qtls.csv", index=False)

qtls_g = qtls_image >> \
            group_trait(property_trait='trait', name='group_category', groups=image_groups) >> \
            group_by(f.group_category,f.chr) >> \
            group_modify(merge_qtl) >> \
            ungroup() >> \
            filter_qtl()
#            filter_qtl(f.count >= 2, f.mean_marker_variance > 10.0)

qtls_g.to_csv("../../Data/publication/tables/diazGarcia2018ImagePhenotyping.supplemental.grouped.consensus.qtls.csv", index=False)

#Merge Schlautman 2015 qtl results with my CNJ0x results
mvariance_min = 5
schlautman_groups=[["MFM","UMFM","UTBM","UBW","UBM","UBL"],["TY","SFY","PAC"]] #Group QTL by traits that have similar meaning/correlation

qtls_schlautman = tibble(pd.read_csv("../../Data/publication/tables/schlautman2015.qtls.csv")) >> \
        rename(interval_left=f.lod_2_min,
               interval_right=f.lod_2_max,
               qtl_lod=f.lod_2) >> \
        filter(((f.interval_right - f.interval_left) <= interval_max) & (f.marker_variance >= mvariance_min)) >> \
        mutate(population="cnj02",
                model_variance=np.nan,
               method="CIM") >> \
        mutate(study="schlautman") >> \
        convert_ssr2composite(ssr2c_class=ssr2c) >> \
        select(columns_highlighted)

nis = nie #Previous value
nie = nis + qtls_schlautman.shape[0]
qtls_schlautman = qtls_schlautman >> \
                mutate(id=range(nis,nie))
qtls_schlautman.index = range(0,len(qtls_schlautman))
qtls_schlautman.to_csv("../../Data/publication/tables/schlautman2015_cnj0x.raw.consensus.csv", index=False)

qtls_merged = qtls_schlautman >> \
                bind_rows(qtls_cnj02) >> \
                bind_rows(qtls_cnj04) >> \
                select(columns_highlighted_i) >> \
                group_trait(property_trait='trait', name='group_category', groups=schlautman_groups) >> \
                mutate(count = 1) >> \
                group_by(f.group_category,f.chr) >> \
                group_modify(merge_qtl) >> \
                ungroup() >> \
                filter_qtl()
#                filter_qtl(f.count >= 6, f.mean_marker_variance > 10.0,f.model_count >= 3)

qtls_merged.to_csv("../../Data/publication/tables/schlautman2015_cnj0x.grouped.consensus.csv", index=False)


#Merge Diaz-Garcia MassivePhenotyping results with CNJ0x results

cnj0x_massive_groups=[["MFM","UMFM","UTBM","UBW","UBM","UBL"],["TY","SFY","PAC"],["Tacy","BCOLOR","BCOLORVAR","TACY_DIFF","TACY_OCT","TACY_SEP"]] #Group QTL by traits that have similar meaning/correlation

qtls_merged = qtls_massive >> \
                bind_rows(qtls_cnj02) >> \
                bind_rows(qtls_cnj04) >> \
                select(columns_highlighted_i) >> \
                group_trait(property_trait='trait', name='group_category', groups=cnj0x_massive_groups) >> \
                mutate(count = 1) >> \
                group_by(f.group_category,f.chr) >> \
                group_modify(merge_qtl) >> \
                ungroup() >> \
                filter_qtl()
#                filter_qtl(f.count >= 6, f.mean_marker_variance > 10.0, f.model_count >= 3)

qtls_merged.to_csv("../../Data/publication/tables/diazGarciaMassivePhenotyping2018_cnj0x.grouped.consensus.csv", index=False)

#Merge Diaz-Garcia ImagePhenotyping results with CNJ0x results
cnj0x_image_groups=[["MFM","UMFM","UTBM","UBW","UBM","UBL","BL","BW","BA"],["TY","SFY","PAC"],["UKEC","ULvW","LvW","EC"]] #Group QTL by traits that have similar meaning/correlation

qtls_merged = qtls_image >> \
                bind_rows(qtls_cnj02) >> \
                bind_rows(qtls_cnj04) >> \
                select(columns_highlighted_i) >> \
                group_trait(property_trait='trait', name='group_category', groups=cnj0x_image_groups) >> \
                mutate(count = 1) >> \
                group_by(f.group_category,f.chr) >> \
                group_modify(merge_qtl) >> \
                ungroup() >> \
                filter_qtl()
#                filter_qtl(f.count >= 6, f.mean_marker_variance > 10.0, f.model_count >= 3)

qtls_merged.to_csv("../../Data/publication/tables/diazGarciaImagePhenotyping2018_cnj0x.grouped.consensus.csv", index=False)

#Merge Schlautman 2015 qtl results and Diaz-Garcia ImagePhenotyping results with CNJ0x results

qtls_merged = qtls_schlautman >> \
                bind_rows(qtls_cnj02) >> \
                bind_rows(qtls_cnj04) >> \
                bind_rows(qtls_image) >> \
                select(columns_highlighted_i) >> \
                group_trait(property_trait='trait', name='group_category', groups=cnj0x_image_groups) >> \
                mutate(count = 1) >> \
                group_by(f.group_category,f.chr) >> \
                group_modify(merge_qtl) >> \
                ungroup() >> \
                filter_qtl()
#                filter_qtl(f.count >= 6, f.mean_marker_variance > 6.0, f.model_count >= 3)

qtls_merged.to_csv("../../Data/publication/tables/schlautman2015_diazGarciaImagePhenotyping2018_cnj0x.grouped.consensus.csv", index=False)

#Dump all ungrouped QTL to a csv file so that IDs can be traced
qtls_raw = qtls_schlautman >> \
                bind_rows(qtls_cnj02) >> \
                bind_rows(qtls_cnj04) >> \
                bind_rows(qtls_image) >> \
                bind_rows(qtls_massive)

qtls_raw.to_csv("../../Data/publication/tables/all.raw.consensus.csv", index=False)
