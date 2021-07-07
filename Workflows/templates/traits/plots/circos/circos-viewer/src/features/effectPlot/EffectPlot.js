import React, { useEffect, useState } from 'react';
import Plot from 'react-plotly.js';
import { useSelector, shallowEqual } from 'react-redux';
import { selectEffectPlotData } from './effectPlotSlice';
import styles from './EffectPlot.module.css';

import * as d3 from 'd3';

const EffectPlot = () => {
    const plotData                              = useSelector(selectEffectPlotData, shallowEqual);
    const [currentPlotData, setCurrentPlotData] = useState({d: [], l: {width: 480, height: 320, title: 'Effect Plot'}});

    useEffect(() => {
        const transformedData = plotData.data.map( (d,i) => {
            const color = d3.hsl(plotData.color);
            const color_hsl = color.copy({h: color.h + (i*25)}); /* Adjust hue slightly to indicate different genotypes */
            const outliercolor_hsl = color_hsl.copy({s: color_hsl.s * 0.75}); /* Drop the saturation to indicate outlier */
            const outliercolor_rgba = outliercolor_hsl.rgb().copy({opacity: 0.6});
            return( {
                x: d.blups.map( b => b.blup ),
                text: d.blups.map( b => ("g"+b.id_i) ),
                type: 'box',
                name: d.genotype,
                jitter: 0.5,
                fillcolor: color_hsl,
                hoverinfo: "text+x",
                marker: {
                    color: d3.color("black"),
                    opacity: 0.5,
                    outliercolor: outliercolor_rgba,
                    size: 3
                },
                line: {
                    color: "black",
                    width: 1
                },
                boxpoints: 'all',
                pointpos: 0
            } );
        } );
        const layout={width: 480, height: 320, title: '<b>Effect Plot</b> ' + plotData.model + ' ' + plotData.trait + '<br>lg' + plotData.chr + '@' + Math.round(plotData.position/1000,1) + 'cM'};
        setCurrentPlotData({d: transformedData, l: layout});
    }, [plotData]);

    return(
        <Plot
            data={currentPlotData.d}
            layout={currentPlotData.l}
        />
    );
}; 

export default EffectPlot;
