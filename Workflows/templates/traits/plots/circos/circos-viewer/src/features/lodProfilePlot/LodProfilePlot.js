import React, { useEffect, useState } from 'react';
import Plot from 'react-plotly.js';
import { useSelector, shallowEqual } from 'react-redux';
import { selectLodProfilePlotData } from './lodProfilePlotSlice';
import styles from './LodProfilePlot.module.css';

import * as d3 from 'd3';

const LodProfilePlot = () => {
    const plotData                              = useSelector(selectLodProfilePlotData, shallowEqual);
    const [currentPlotData, setCurrentPlotData] = useState({d: [], l: {width: 480, height: 320, title: '<b>LOD Profile Plot</b>'}});

    useEffect(() => {
        const data = {
            x: plotData.data.map( d => d.position ),
            y: plotData.data.map( d => d.lod ),
            type: 'scatter',
            mode: 'lines+markers',
            marker: {
                color: d3.color("black"),
                opacity: 0.5,
                size: 3
            },
            line: {
                shape: 'linear',
                color: plotData.color,
                width: 1
            },
            fill: "tozeroy",
            fillColor: d3.color(plotData.color).rgb().copy({opacity: 0.75}),
        };
        const layout={...currentPlotData.l, title: '<b>LOD Profile Plot</b>: ' + plotData.model + ' ' + plotData.trait + '<br>LG' + plotData.chr};
        setCurrentPlotData({d: [data], l: layout});
    }, [plotData]);

    return(
        <Plot
            data={currentPlotData.d}
            layout={currentPlotData.l}
        />
    );
}; 

export default LodProfilePlot;
