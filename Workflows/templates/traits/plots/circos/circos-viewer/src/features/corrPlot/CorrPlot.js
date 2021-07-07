import React, { useEffect, useState } from 'react';
import Plot from 'react-plotly.js';
import { useSelector, shallowEqual } from 'react-redux';
import { selectCorrPlotTraits } from './corrPlotSlice';
import {selectList} from '../viewController/viewControllerSlice';
import styles from './CorrPlot.module.css';

import * as d3 from 'd3';
import DataFrame from 'dataframe-js';

const CorrPlot = () => {
    const [init, setInit]           = useState(false);
    const [plotData, setPlotData]   = useState(undefined);
    const [traitmap, setTraitmap]   = useState(undefined);
    const traits                    = useSelector(selectList('traits'));
    const [currentPlotData, setCurrentPlotData] = useState({d: [], l: {width: 640, height: 320, title: '<b>Correlation Heatmap</b>'}});

    useEffect(() => {
        d3.csv('configs/corrplot.cnj02.csv')
          .then( d => {
              setPlotData(new DataFrame(d));
          });
        d3.csv('configs/model-traits.cfg.csv')
          .then( d => {
              var configs_df = new DataFrame(d);
              var traitmap_obj = {};
              configs_df.reduce( (p,n) => (traitmap_obj[n.get('trait')] = {label: n.get('label'), short: n.get('label_short')}), {});
              setTraitmap(traitmap_obj);
          });
    }, [init]);
    if( init == false ) {
        setInit(true);
    }

    useEffect(() => {
        if( (plotData != undefined) && (traitmap != undefined) && (traits.length > 0) ) {
            const traits_short = traits.filter( t => t.enabled ).map( t => traitmap[t.id].short );
            const traits_long = traits.filter( t => t.enabled ).map( t => t.text );
            const columns = ['',...traits_short]; //The '' column represents the row names
            var dataSubset = plotData.select(...columns).filter( r => traits_short.includes(r.get('')) );
            var dataSubsetA = dataSubset.drop('').toArray();
            for( var i = 0; i < dataSubsetA.length; i++ ) {
                var row = dataSubsetA[i];
                let rname = columns[i+1];
                let ridx  = columns.indexOf(rname)-1;
                for( var j = ridx; j < row.length; j++ ) {
                    row[j] = null;
                }
            };
            const transformedData = {
                x: traits_long,
                y: traits_long.reverse(),
                z: dataSubsetA.reverse(),
                type: 'heatmap',
                hoverongaps: false
            };
            setCurrentPlotData({...currentPlotData, d: [transformedData]});
        }
    }, [plotData, traits, traitmap]);

    return(
        <Plot
            data={currentPlotData.d}
            layout={currentPlotData.l}
        />
    );
}; 

export default CorrPlot;
