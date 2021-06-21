import React, {useState, useEffect} from 'react';
import { useDispatch, useSelector } from 'react-redux';
import {selectQTLConsensus,
        selectQTLMethod,
        selectList,
        selectDisplayTrackLabels,
        setList
} from '../viewController/viewControllerSlice';

import * as d3 from 'd3';
import DataFrame from 'dataframe-js';

import Circos from 'circos';

var qtl_bubble_opacity            = 0.75;
var qtl_bubble_highlight_opacity  = 1.0;
var qtl_bubble_scale_factor       = 15;
var circosScatter = undefined;
var circos_trait2traitname = {};
var circos_trait2traitshort = {};
var circos_model2modelname = {};
var gkaryotypes;
var glayout;
var gtrait_files;
var glods;
var gqtls;
var gconfigs;
var gcircos_trait_config;
var gscatter_defaults;
var gstack_defaults;
var gline_defaults;
var blipColors = d3.scaleSequential(d3.interpolateTurbo);
var axesColors = d3.scaleSequential(d3.interpolateGreys);

const seq = (start, stop, step) => Array.from({ length: (stop - start) / step + 1}, (_, i) => start + (i * step));

var change_scan_type = function(selected_type) {
    //Hide everything
    //Need to hide any old state for the LOD profile lines.
    d3.selectAll("*[class^='lines-']")
        .attr('visibility', 'hidden');
    //Hide all current scatters-
    d3.selectAll("*[class^='scatters-']")
        .attr('opacity','0')
        .attr('visibility','hidden');
    d3.selectAll("*[class^='stacks-']")
        .attr('opacity','0')
        .attr('visibility','hidden');

    //Show what we care about
    d3.selectAll("*[class^='scatters-"+selected_type+"--']")
        .attr('opacity','1')
        .attr('visibility','visible');
    d3.selectAll("*[class^='stacks-"+selected_type+"--']")
        .attr('opacity','1')
        .attr('visibility','visible');
}

var get_current_config = function(config,callback) {
    setTimeout(function() {
        callback(null, config);
    }, 250);
}

var show_lod_profs = function(datum, index, nodes, event) {
    var mclass = 'lines-'+datum.method+'--'+datum.model+'--'+datum.trait;
    //Hide other lodprofiles that are displaying in the current trait
    d3.selectAll("*[class^='lines-'][class$='"+datum.trait+"']")
        .attr('visibility', 'hidden');
    d3.selectAll("*[class='"+mclass+"']")
        .attr('visibility', 'show');
}

var hide_lod_profs = function(datum, index, nodes, event) {
    var mclass = 'lines-'+datum[0].method+'--'+datum[0].model+'--'+datum[0].trait;
    d3.selectAll("*[class='"+mclass+"']")
        .attr('visibility', 'hidden');
}

var highlight_markers = function(datum, index, nodes, event) {
    var marker = datum.nearest_marker;
    d3.selectAll("*[marker='"+marker+"']")
        .attr('opacity', qtl_bubble_highlight_opacity);
}

var unhighlight_markers = function(datum, index, nodes, event) {
    var marker = datum.nearest_marker;
    d3.selectAll("*[marker='"+marker+"']")
        .attr('opacity', qtl_bubble_opacity);
}

const inject_scatters = (linkage_groups, qtls, traits, trait_idx, models, config, method, consensus) => {
  var trait = traits[trait_idx];
  var scatter_trait_data = qtls.filter( d => (d.trait === trait) )
                                .filter( d => (models.includes(d.model)) )
                                .filter( d => ((d.method === method) || (d.method === "fakeqtl")) )
                                .map( d => ({   ...d,
                                                block_id: "vm"+d.chr,
                                                color: blipColors(((trait_idx * models.length) + models.indexOf(d.model))/(traits.length * models.length)),
                                                position: +(consensus === "consensus" ? d.position_consensus : d.position) * 1000,
                                                value: models.indexOf(d.model),
                                                size: Math.round(d.marker_variance * qtl_bubble_scale_factor),
                                                lod: +d.qtl_lod,
                                                pval: +d.qtl_pvalue,
                                                class: d.class
                                          }) )
                                .filter( d => (linkage_groups.includes(d.block_id)) );
  var trait_class = "scatters-" + method + "-" + consensus + "--" + trait;
  circosScatter.scatter(trait_class, scatter_trait_data, config);
};

const inject_stacks = (linkage_groups, qtls, traits, trait_idx, models, config, method, consensus) => {
    var trait = traits[trait_idx];
    models.forEach( (model,i) => {
        //scanone normal
        var stack_trait_data = qtls.filter( d => (d.trait === trait) )
                                    .filter( d => (d.model === model) )
                                    .filter( d => ((d.method === "scanone") || (d.method === "fakeqtl")) )
                                    .map( d => ({ ...d,
                                                    block_id: "vm" + d.chr,
                                                    color: blipColors(((trait_idx * models.length) + models.indexOf(d.model))/(traits.length * models.length)),
                                                    start: (+(consensus === "consensus" ? d.position_consensus : d.position) - (+d.interval/2))*1000,
                                                    end: (+(consensus === "consensus" ? d.position_consensus : d.position) + (+d.interval/2))*1000
                                                }))
                                    .filter( d => (linkage_groups.includes(d.block_id)) );
        if( stack_trait_data.length > 0 ) {
            var trait_class = "stacks-"+method+"-"+consensus+"--"+model+"--"+trait;
            circosScatter.stack(trait_class, stack_trait_data, config[i]);
        } else {
            console.log("Non-zero stack_trait_data length.")
        }
    });
};

const inject_lods = (linkage_groups, trait, models, config, method) => {
    var lconfig = {...config}; //Make a copy
    models.forEach( (model,i) => {
        if( model in glods[trait] ) {
            var lods_trait_model = glods[trait][model];
            if( method in lods_trait_model ) {
                var lods_method = lods_trait_model[method];
                linkage_groups.forEach( (lg,j) => {
                    if( lg in lods_method ) {
                        var lod_trait_data = [];
                        var lods = lods_method[lg];
                        lods.forEach( d => {
                            lod_trait_data.push( {
                                method: method,
                                block_id: lg,
                                position: +d.position*1000,
                                value: +d.lod,
                                model: model,
                                trait: trait
                            });
                        } );
                        //Edit the Lconfig  max lod values, axes spacing, etc.
                        var max_lod = lod_trait_data.reduce( (a,b) => (a.value > b.value ? a : b) ).value;
                        lconfig.max = max_lod + max_lod/10; // Add 10% to top
                        lconfig.backgrounds[0].end = lconfig.max;
                        lconfig.axes[0].spacing = lconfig.max/(models.length);
                        var trait_class = "lines-"+method+"--"+model+"--"+trait;
                        circosScatter.line(trait_class, lod_trait_data, lconfig);
                    }
                });
            }
        }
    });
}

const inject_data = (linkage_groups, qtls, traits, trait_idx, models, scatter_config, stack_config, line_config) => {
    //Scatterplots (QTL bubble tracks)
    inject_scatters(linkage_groups, qtls, traits, trait_idx, models, scatter_config, "scanone", "normal");
    inject_scatters(linkage_groups, qtls, traits, trait_idx, models, scatter_config, "scanone", "consensus");
    inject_scatters(linkage_groups, qtls, traits, trait_idx, models, scatter_config, "stepwiseqtl", "normal");
    inject_scatters(linkage_groups, qtls, traits, trait_idx, models, scatter_config, "stepwiseqtl", "consensus");

    //Stacks (1.5 LOD interval tracks)
    inject_stacks(linkage_groups, qtls, traits, trait_idx, models, stack_config, "scanone", "normal");
    inject_stacks(linkage_groups, qtls, traits, trait_idx, models, stack_config, "scanone", "consensus");
    inject_stacks(linkage_groups, qtls, traits, trait_idx, models, stack_config, "stepwiseqtl", "normal");
    inject_stacks(linkage_groups, qtls, traits, trait_idx, models, stack_config, "stepwiseqtl", "consensus");

    //LOD profiles
    inject_lods(linkage_groups, traits[trait_idx], models, line_config, "scanone");
    inject_lods(linkage_groups, traits[trait_idx], models, line_config, "stepwiseqtl");
}

const gen_scatter_color = d => (d.color);

const gen_scatter_size = d => (d.size);

const gen_scatter_tooltip = d => {
    var variance = Math.round(d.marker_variance, 1);
    var position = Math.round(d.position/1000.0, 1);
    return("<b>"+circos_trait2traitname[d.trait]+"</b><br />"+circos_model2modelname[d.model]+"<br />Pos: "+position.toString()+"cM<br />Var: "+variance.toString()+"%<br />"+d.nearest_marker);
}

const gen_stack_tooltip = function(d) {
    var start = Math.round(d.start/1000.0, 1);
    var end = Math.round(d.end/1000.0, 1);
    var range = Math.round((d.end-d.start)/1000.0, 1);
    return("<b>"+circos_trait2traitname[d.trait]+"</b><br />"+circos_model2modelname[d.model]+"<br />Start-End: "+start.toString()+"-"+end.toString()+"cM<br />Range: "+range.toString()+"cM");
}

const drawCircos = (linkage_groups, qtls, traits, models, track_configs) => {
  var relevant_karyotypes = gkaryotypes.filter( k => (linkage_groups.includes(k.id)) );
  circosScatter
    .layout(
        relevant_karyotypes,
        glayout
    );
  for( var i = 0; i < traits.length; i++ ) {
    inject_data(linkage_groups, qtls, traits, i, models, track_configs.scatters[i], track_configs.stacks[i], track_configs.lods[i]);
  }
  circosScatter.render();
  d3.selectAll(".point")
    .on('mouseover.highlight', highlight_markers)
    .on('mouseout.highlight', unhighlight_markers);
  d3.selectAll(".tile")
    .on('mouseover.highlight', highlight_markers)
    .on('mouseout.highlight', unhighlight_markers);
  d3.selectAll("g[class='line']")
      .on('click', hide_lod_profs);
  d3.selectAll("*[class^='lines-']")
      .attr('opacity','0.8')
      .attr('visibility','hidden');
  d3.selectAll("*[class^='scatters-']")
      .attr('opacity','0')
      .attr('visibility','hidden');
  d3.selectAll("*[class^='stacks-']")
      .attr('opacity','0')
      .attr('visibility','hidden');
  //Find the arc of the centermost scatterplot and build a square.
          d3.selectAll("*[class^='scatter']")
          .selectAll(".block")
          .select("path.background")
          .attr("d");
                                  
  //Update all the scatter plots with a 'marker' attribute to select all identical nearest markers
  var points = d3.selectAll(".point")
                 .datum( function(d) {
                     return d;
                 })
                 .attr( "marker", function(d) {
                      return d.nearest_marker;
                 });
  //Update all stack intervals with a 'marker' attribute to select all identical nearest interval stacks
  var stacks = d3.selectAll(".tile")
                 .datum( function(d) {
                     return d;
                 })
                 .attr( "marker", function(d) {
                      return d.nearest_marker;
                 });
  d3.selectAll("*[marker^='marker--fake--trait']").attr('visibility', 'hidden');
}

const globalTrackGenerator = (ntraits, intertrack_distance=0.3, track_start=0.3, track_end=1.0) => {
  let gband   = d3.scaleBand()
                  .domain(seq(0,ntraits,1))
                  .range([track_start, track_end])
                  .paddingInner(intertrack_distance);
  return(gband);
};

const stackTrackGenerator = (gband, trait_index, nmodels, intertrack_distance=0.01, stack_proportion=0.3) => {
  const trait_start = gband(trait_index) + gband.bandwidth();
  const trait_end = gband(trait_index+1);
  let sband = d3.scaleBand()
                .domain(seq(0, nmodels, 1))
                .range([trait_start, trait_end])
                .padding(intertrack_distance);
  return(sband);
};

const reloadTrackConfigs = (linkage_groups, traits, models) => {
  //Loop through the traits and generate scatter configs
  var ntraits = traits.length;
  var scatter_configs_json = Array(ntraits);
  var stack_configs_json = Array(ntraits);
  var line_configs_json = Array(ntraits);
  var nmodels = models.length;
  const trackBands = globalTrackGenerator(ntraits);
  for( let i = 0; i < ntraits; i++ ) {
      var trait = traits[i];
      //Edit the scatter layout
      var scatter_config_json = JSON.parse(JSON.stringify(gscatter_defaults)); //Deep copy
      var bi = trackBands(i);
      scatter_config_json.innerRadius     = bi;
      scatter_config_json.outerRadius     = bi + trackBands.bandwidth();
      scatter_config_json.max             = nmodels;
      scatter_config_json.color           = gen_scatter_color;
      scatter_config_json.size            = gen_scatter_size;
      scatter_config_json.tooltipContent  = gen_scatter_tooltip;
      scatter_config_json.trackLabelConf.label = circos_trait2traitname[trait];

      //Generate the axes
      var axis_template                    = scatter_config_json.axes[0];
      scatter_config_json.axes         = Array(nmodels);
      for( let j = 0; j < nmodels; j++ ) {
          scatter_config_json.axes[j]                       = JSON.parse(JSON.stringify(axis_template));
          scatter_config_json.axes[j].position              = j+1;
          scatter_config_json.axes[j].color                 = axesColors((j+1)/nmodels);
          scatter_config_json.axes[j].axisLabelConf.label   = models[j];
      }
      scatter_config_json.backgrounds[0].start    = 0;
      scatter_config_json.backgrounds[0].end      = nmodels;
      scatter_config_json.backgrounds[0].color    = d3.color("gray45");
      scatter_config_json.backgrounds[0].opacity  = 1;
      scatter_configs_json[i] = scatter_config_json;
      //Edit the stack layout
      var stackBands = stackTrackGenerator(trackBands, i, nmodels);
      var stack_configs_inner_json = Array(nmodels);
      for( let j = 0; j < nmodels; j++ ) {
          var stack_config_json = JSON.parse(JSON.stringify(gstack_defaults));
          stack_config_json.innerRadius  = stackBands(j);
          stack_config_json.outerRadius  = stackBands(j) + stackBands.bandwidth();
          stack_config_json.backgrounds[0].start = 0;
          stack_config_json.backgrounds[0].color = d3.color("gray22");
          stack_config_json.backgrounds[0].opacity = 1;
          stack_config_json.color = gen_scatter_color;
          stack_config_json.tooltipContent = gen_stack_tooltip;
          stack_configs_inner_json[j] = stack_config_json;
      }
      stack_configs_json[i] = stack_configs_inner_json;
  
      //LOD Line Configurations
      var line_config_json = JSON.parse(JSON.stringify(gline_defaults));
      line_config_json.innerRadius  = bi;
      line_config_json.outerRadius  = bi + trackBands.bandwidth();
      line_config_json.color = "white";
      line_config_json.backgrounds[0].color = d3.color("gray45");
      line_configs_json[i] = line_config_json;
  }
  return( { scatters: scatter_configs_json,
            stacks: stack_configs_json,
            lods: line_configs_json } );
};

const generateNameMaps = (tconfig) => {
  for( var i = 0; i < tconfig.length; i++ ) {
    var trait = tconfig[i].trait;
    var model = tconfig[i].model;
    if(!(trait in circos_trait2traitname)) {
        circos_trait2traitname[trait] = tconfig[i].label;
    }
    if(!(trait in circos_trait2traitshort)) {
        circos_trait2traitshort[trait] = tconfig[i].label_short;
    }
    if(!(model in circos_model2modelname)) {
        circos_model2modelname[model] = tconfig[i].model_label;
    }
  }
}

const refresh = (linkage_groups, qtls, traits, models) => {
  gconfigs = reloadTrackConfigs(linkage_groups,traits,models);
  drawCircos(linkage_groups, qtls, traits, models, gconfigs);
}

const inject_fake_qtls = (layout, linkage_groups, traits, models, qtls) => {
    const fake_qtls = [];
    var trackLabelBlockId = "";
    if( 'trackLabelBlockId' in layout ) {
        trackLabelBlockId = layout.trackLabelBlockId;
    }
    linkage_groups.forEach( lg => {
        if( lg !== trackLabelBlockId ) { //Don't inject data on the label track -- just want the text to render
            traits.forEach( t => {
                models.forEach( m => {
                    fake_qtls.push( {
                        method: "fakeqtl",
                        chr: lg.replace( /(vm)(\d+)$/, '$2' ),
                        trait: t,
                        model:  m,
                        position: -1,
                        position_consensus: -1,
                        chr2: undefined,
                        position2: undefined,
                        nearest_marker: "marker--fake--trait--" + t,
                        qtl_lod: 0,
                        qtl_pvalue: 0,
                        marker_variance: 0,
                        model_variance: 0,
                        interval: 0,
                        class: "" 
                    } );
                } );
            } );
        }
    } );
    fake_qtls.push(...qtls);
    return(fake_qtls);
}


const generate_new_circos = (container="#scatterChart", width="1024", height="1024", opacity=0.75, highlight_opacity=1.0, bubble_scale=15) => {
    if( circosScatter !== undefined ) {
        circosScatter.detach();
        delete window.circosScatter;
    }
    circosScatter = Circos({
        container: container,
        width: width,
        height: height 
    });

    qtl_bubble_opacity            = opacity;
    qtl_bubble_highlight_opacity  = highlight_opacity;
    qtl_bubble_scale_factor       = bubble_scale;
}

export default function View() {
    const dispatch                    = useDispatch(); 
    const qtlConsensus                = useSelector(selectQTLConsensus);
    const qtlMethod                   = useSelector(selectQTLMethod);
    const displayTrackLabels          = useSelector(selectDisplayTrackLabels);
    const linkageGroups               = useSelector(selectList('linkageGroups'));
    const models                      = useSelector(selectList('models'));
    const traits                      = useSelector(selectList('traits'));
    const [init, setInit]             = useState(false);
    const [postInit, setPostInit]     = useState(false);

    const generateMethodConsensus = () => {
        var consensus = qtlConsensus ? "consensus" : "normal";
        return(qtlMethod + '-' + consensus);
    };

    const loadStaticConfigurations = (data) => {
        //Copy the loaded static configuration to globals
        gkaryotypes           = data[0];
        glayout               = data[1];
        gcircos_trait_config  = data[2];
        gqtls                 = data[3];
        glods                 = data[4];
        gscatter_defaults     = data[5];
        gstack_defaults       = data[6];
        gline_defaults        = data[7];

        /* By default, enable everything */
        var karyotypes_df  = new DataFrame(gkaryotypes); 
        var config_df      = new DataFrame(gcircos_trait_config);
        var linkage_groups = karyotypes_df.select('id', 'label');
        dispatch(setList('linkageGroups')(linkage_groups.toCollection().map(e => ({id: e.id, text: e.label, enabled: true}))));
        var models         = config_df.unique('model').toArray('model');
        var model_names    = config_df.unique('model_label').toArray('model_label');
        dispatch(setList('models')(models.map((m,i) => ({id: m, text: model_names[i], enabled: true}))));
        var traits         = config_df.unique('trait').toArray('trait');
        var trait_names    = config_df.unique('label').toArray('label');
        dispatch(setList('traits')(traits.map((t,i) => ({id: t, text: trait_names[i], enabled: true}))));

        gqtls = inject_fake_qtls(glayout, linkage_groups.select('id').toArray('id'), traits, models, gqtls);

        generateNameMaps(gcircos_trait_config);
    };

    useEffect(() => {
        /* This is where I should load static configuration */
        /*
        */
        Promise.all([   d3.json('configs/circos/karyotype.json'),
                        d3.json('configs/circos/layout.default.json'),
                        d3.csv('configs/model-traits.cfg.csv'),
                        d3.csv('configs/circos/qtl_collated.consensus.csv'),
                        d3.json('configs/circos/lod_profiles.json'),
                        d3.json('configs/circos/scatter.default.json'),
                        d3.json('configs/circos/stack.default.json'),
                        d3.json('configs/circos/line.default.json')])
                .then(data => {
                    loadStaticConfigurations(data);
                    setPostInit(true);
                })
                .catch(error => {
                    console.error(error.message)
                });
        return;
    }, [init]);

    useEffect(() => {
        change_scan_type(generateMethodConsensus());
    }, [qtlConsensus, qtlMethod]);

    useEffect(() => {
        if( postInit === true ) {
            var labelSector = glayout['trackLabelBlockId']
            var linkage_groups    = linkageGroups
                                        .filter(k => k.enabled)
                                        .map(k => k.id);
            var traits_enabled    = traits
                                        .filter(t => t.enabled)
                                        .map(t => t.id);
            var models_enabled    = models
                                        .filter(m => m.enabled)
                                        .map(m => m.id);
            if( displayTrackLabels ) {
                gkaryotypes.push(
                {
                    "id": labelSector,
                    "label": "", //Must be blank
                    "color": "rgb(255,255,255)",
                    "len": 230000
                }); //Inject a sector to show labels
                linkage_groups.push(labelSector);
            } else {
                gkaryotypes = gkaryotypes.filter(k => (k.id !== labelSector));
            }
            generate_new_circos();
            refresh(linkage_groups, gqtls, traits_enabled, models_enabled);
            change_scan_type(generateMethodConsensus());
        }
    }, [postInit, displayTrackLabels, linkageGroups, traits, models]);


    if( init === false ) {
        setInit(true);
    }

    return (
        <div id='scatterChart'></div>
    );
}
