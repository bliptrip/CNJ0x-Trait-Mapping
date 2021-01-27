var qtl_bubble_opacity            = 0.75;
var qtl_bubble_highlight_opacity  = 1.0;
var qtl_bubble_scale_factor       = 15;
var gfolder_prefix                = "./";
var circosScatter;
var circos_trait2traitname = {};

var change_scan_type = function(selected_type) {
    //Hide everything
    //Need to hide any old state for the LOD profile lines.
    d3.selectAll("*[class^='lines-']")
        .attr('visibility', 'hidden');
    //Hide all current scatter-
    d3.selectAll("*[class^='scatter-']")
        .attr('opacity','0')
        .attr('visibility','hidden');
    d3.selectAll("*[class^='stacks-']")
        .attr('opacity','0')
        .attr('visibility','hidden');

    //Show what we care about
    d3.selectAll("*[class^='scatter-"+selected_type+"--']")
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

var inject_scatter = function(error, data)
{
    bin_size = 3; //Specifies the number of indices per 'bin' of data
    for( i = 0; i < (data.length/bin_size); i++ ) {
        const trait_data = data[bin_size*i];
        var   models = [];
        const map = new Map(); //For tracking distinct elements
        for (const d of trait_data) {
                if(!map.has(d.model)){
                        map.set(d.model, true);    // Only add if not already added
                        models.push({
                                        model: d.model,
                                        model_idx: d.model_idx
                                   });
                }
        }
        models = models.sort( (a,b) => {
                                if(a.model_idx > b.model_idx) {
                                        return(1);
                                } else if(a.model_idx < b.model_idx) {
                                        return(-1);
                                } else {
                                        return(0);
                                }
                        });
        model_names = models.map( d => d.model );
        scatter_config     = data[(bin_size*i)+1];
        scatter_config.trackLabelConf.label = circos_trait2traitname[trait_data[0].trait]
        scatter_config.axes.map( (a,i) => {
            var b = {...a};
            b.axisLabelConf.label = model_names[i];
            return(b);
        });
        scatter_config.events = {
            'click.showlodprofs': show_lod_profs
        };
        //scanone results
        scatter_trait_data = trait_data
                                     .filter(function(d) { return((d.method == "scanone") || (d.method == "fakeqtl")); })
                                     .map(function(d) {
                                             return({
                                                     method: d.method,
                                                     block_id: "vm"+d.chr,
                                                     position: +d.position_consensus*1000,
                                                     value: +d.model_idx,  
                                                     model: d.model,
                                                     trait: d.trait,
                                                     color: d.color,
                                                     stroke_color: d.stroke_color,
                                                     model_color: d.model_color,
                                                     nearest_marker: d.nearest_marker,
                                                     variance: d.marker_variance,
                                                     size: Math.round(d.marker_variance * qtl_bubble_scale_factor),
                                                     lod: +d.qtl_lod,
                                                     pval: +d.qtl_pvalue,
                                                     class: d.class
                                             });
                                     });
        trait = "scatter-scanone-consensus--" + scatter_trait_data[0].trait;
        circosScatter.scatter(trait, scatter_trait_data, scatter_config);
        scatter_trait_data = trait_data
                                     .filter(function(d) { return((d.method == "scanone") || (d.method == "fakeqtl")); })
                                     .map(function(d) {
                                             return({
                                                     method: d.method,
                                                     block_id: "vm"+d.chr,
                                                     position: +d.position*1000,
                                                     value: +d.model_idx,  
                                                     model: d.model,
                                                     trait: d.trait,
                                                     color: d.color,
                                                     stroke_color: d.stroke_color,
                                                     model_color: d.model_color,
                                                     nearest_marker: d.nearest_marker,
                                                     variance: d.marker_variance,
                                                     size: Math.round(d.marker_variance * qtl_bubble_scale_factor),
                                                     lod: +d.qtl_lod,
                                                     pval: +d.qtl_pvalue,
                                                     class: d.class
                                             });
                                     });
        trait = "scatter-scanone-normal--" + scatter_trait_data[0].trait;
        circosScatter.scatter(trait, scatter_trait_data, scatter_config);
        //stepwiseqtl results
        scatter_trait_data = trait_data
                                     .filter(function(d) { return((d.method == "stepwiseqtl") || (d.method == "fakeqtl")); })
                                     .map(function(d) {
                                             return({
                                                     method: d.method,
                                                     block_id: "vm"+d.chr,
                                                     position: +d.position_consensus*1000,
                                                     value: +d.model_idx,  
                                                     model: d.model,
                                                     trait: d.trait,
                                                     color: d.color,
                                                     stroke_color: d.stroke_color,
                                                     model_color: d.model_color,
                                                     nearest_marker: d.nearest_marker,
                                                     variance: d.marker_variance,
                                                     size: Math.round(d.marker_variance * qtl_bubble_scale_factor),
                                                     lod: +d.qtl_lod,
                                                     pval: +d.qtl_pvalue,
                                                     class: d.class
                                             });
                                     });
        trait = "scatter-stepwiseqtl-consensus--" + scatter_trait_data[0].trait;
        circosScatter.scatter(trait, scatter_trait_data, scatter_config);
        scatter_trait_data = trait_data
                                     .filter(function(d) { return((d.method == "stepwiseqtl") || (d.method == "fakeqtl")); })
                                     .map(function(d) {
                                             return({
                                                     method: d.method,
                                                     block_id: "vm"+d.chr,
                                                     position: +d.position*1000,
                                                     value: +d.model_idx,  
                                                     model: d.model,
                                                     trait: d.trait,
                                                     color: d.color,
                                                     stroke_color: d.stroke_color,
                                                     model_color: d.model_color,
                                                     nearest_marker: d.nearest_marker,
                                                     variance: d.marker_variance,
                                                     size: Math.round(d.marker_variance * qtl_bubble_scale_factor),
                                                     lod: +d.qtl_lod,
                                                     pval: +d.qtl_pvalue,
                                                     class: d.class
                                             });
                                     });
        trait = "scatter-stepwiseqtl-normal--" + scatter_trait_data[0].trait;
        circosScatter.scatter(trait, scatter_trait_data, scatter_config);
        stack_configs     = data[(bin_size*i)+2];
        stack_configs     = stack_configs.map( (s,i) => {
            var sc = {...s};
            sc.trackLabelConf.label = model_names[i];
            return(sc);
        });
        for( j = 0; j < stack_configs.length; j++ ) {
            //scanone normal
            stack_trait_data = trait_data
                .filter(function(d) 
                { 
                    var match = (+d.model_idx == (j+1)) && ((d.method == "scanone") || (d.method == "fakeqtl"));
                    return(match);
                })
                .map(function(d) {
                    return({
                        block_id: "vm"+d.chr,
                        start: (+d.position - (+d.interval/2))*1000,
                        end: (+d.position + (+d.interval/2))*1000,
                        model: d.model,
                        trait: d.trait,
                        color: d.color,
                        stroke_color: d.stroke_color,
                        model_color: d.model_color,
                        nearest_marker: d.nearest_marker,
                        class: d.class
                    });
                });
            if( stack_trait_data.length > 0 ) {
                trait = "stacks-scanone-normal--"+stack_trait_data[0].model+"--" + stack_trait_data[0].trait;
                circosScatter.stack(trait, stack_trait_data, stack_configs[j]);
            } else {
                console.log("Non-zero stack_trait_data length.")
            }
            //scanone consensus 
            stack_trait_data = trait_data
                .filter(function(d) 
                { 
                    var match = (+d.model_idx == (j+1)) && ((d.method == "scanone") || (d.method == "fakeqtl"));
                    return(match);
                })
                .map(function(d) {
                    return({
                        block_id: "vm"+d.chr,
                        start: (+d.position_consensus - (+d.interval/2))*1000,
                        end: (+d.position_consensus + (+d.interval/2))*1000,
                        model: d.model,
                        trait: d.trait,
                        color: d.color,
                        stroke_color: d.stroke_color,
                        model_color: d.model_color,
                        nearest_marker: d.nearest_marker,
                        mposition: +d.position_consensus,
                        class: d.class
                    });
                });
            if( stack_trait_data.length > 0 ) {
                trait = "stacks-scanone-consensus--" + stack_trait_data[0].model + "--" + stack_trait_data[0].trait;
                circosScatter.stack(trait, stack_trait_data, stack_configs[j]);
            } else {
                console.log("Non-zero stack_trait_data length.")
            }
            //stepwiseqtl normal
            stack_trait_data = trait_data
                .filter(function(d) 
                { 
                    var match = (+d.model_idx == (j+1)) && ((d.method == "stepwiseqtl") || (d.method == "fakeqtl"));
                    return(match);
                })
                .map(function(d) {
                    return({
                        block_id: "vm"+d.chr,
                        start: (+d.position - (+d.interval/2))*1000,
                        end: (+d.position + (+d.interval/2))*1000,
                        model: d.model,
                        trait: d.trait,
                        color: d.color,
                        stroke_color: d.stroke_color,
                        model_color: d.model_color,
                        nearest_marker: d.nearest_marker,
                        class: d.class
                    });
                });
            if( stack_trait_data.length > 0 ) {
                trait = "stacks-stepwiseqtl-normal--"+stack_trait_data[0].model+"--"+ stack_trait_data[0].trait;
                circosScatter.stack(trait, stack_trait_data, stack_configs[j]);
            } else {
                console.log("Non-zero stack_trait_data length.")
            }
            //stepwiseqtl consensus 
            stack_trait_data = trait_data
                .filter(function(d) 
                { 
                    var match = (+d.model_idx == (j+1)) && ((d.method == "stepwiseqtl") || (d.method == "fakeqtl"));
                    return(match);
                })
                .map(function(d) {
                    return({
                        block_id: "vm"+d.chr,
                        start: (+d.position_consensus - (+d.interval/2))*1000,
                        end: (+d.position_consensus + (+d.interval/2))*1000,
                        model: d.model,
                        trait: d.trait,
                        color: d.color,
                        stroke_color: d.stroke_color,
                        model_color: d.model_color,
                        nearest_marker: d.nearest_marker,
                        mposition: +d.position_consensus,
                        class: d.class
                    });
                });
            if( stack_trait_data.length > 0 ) {
                trait = "stacks-stepwiseqtl-consensus--"+stack_trait_data[0].model+"--" + stack_trait_data[0].trait;
                circosScatter.stack(trait, stack_trait_data, stack_configs[j]);
            } else {
                console.log("Non-zero stack_trait_data length.")
            }
        }
    }
}

var inject_line = function(error, data)
{
    bin_size = 2; //Specifies the number of indices per 'bin' of data
    for( i = 0; i < (data.length/bin_size); i++ ) {
        lod_data     = data[bin_size*i];
        lod_config   = data[(bin_size*i)+1];
        for( var model in lod_data ) {
            //scanone
            lod_trait_data = lod_data[model]
                .filter(function(d) { return(d.method == "scanone"); })
                .map(function(d) {
            //I need to figure out how to separate the different models into different datasets
                    return({
                        method: d.method,
                        block_id: "vm"+d.chr,
                        position: +d.position*1000,
                        value: +d.lod,
                        model: d.model,
                        trait: d.trait
                    });
                });
            //Edit the lod_config per the max lod values, axes spacing, etc.
            max_lod = lod_trait_data.reduce(function(a,b) { 
                return(a.value > b.value ? a : b);
            }).value;
            lod_config.max = max_lod + max_lod/10; // Add 10% to top
            lod_config.backgrounds[0].end = lod_config.max;
            //Have 4 horizontal axes per dataset
            lod_config.axes[0].spacing = lod_config.max/4;
            trait = "lines-scanone--"+model+"--"+lod_trait_data[0].trait;
            circosScatter.line(trait, lod_trait_data, lod_config);
            //stepwiseqtl
            lod_trait_data = lod_data[model]
                .filter(function(d) { return(d.method == "stepwiseqtl"); })
                .map(function(d) {
            //I need to figure out how to separate the different models into different datasets
                    return({
                        method: d.method,
                        block_id: "vm"+d.chr,
                        position: +d.position*1000,
                        value: +d.lod,
                        model: d.model,
                        trait: d.trait
                    });
                });
            //Edit the lod_config per the max lod values, axes spacing, etc.
            max_lod = lod_trait_data.reduce(function(a,b) { 
                return(a.value > b.value ? a : b);
            }).value;
            lod_config.max = max_lod + max_lod/10; // Add 10% to top
            lod_config.backgrounds[0].end = lod_config.max;
            //Have 4 horizontal axes per dataset
            lod_config.axes[0].spacing = lod_config.max/4;
            trait = "lines-stepwiseqtl--"+model+"--"+lod_trait_data[0].trait;
            circosScatter.line(trait, lod_trait_data, lod_config);
        }
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
    d3.selectAll("*[class^='scatter-']")
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

var map_scatter = function(data) {
    return(data);
}

var gen_scatter_color = function(d) {
    return d.color;
}

var gen_scatter_stroke_color = function(d) {
    return d.stroke_color;
}

var gen_scatter_tooltip = function(d) {
        var variance = Math.round(d.variance, 1);
        var position = Math.round(d.position/1000.0, 1);
        return("<b>"+circos_trait2traitname[d.trait]+"</b><br />"+d.model+"<br />Pos: "+position.toString()+"cM<br />Var: "+variance.toString()+"%<br />"+d.nearest_marker);
}

var gen_stack_tooltip = function(d) {
        var start = Math.round(d.start/1000.0, 1);
        var end = Math.round(d.end/1000.0, 1);
        var range = Math.round((d.end-d.start)/1000.0, 1);
        return("<b>"+circos_trait2traitname[d.trait]+"</b><br />"+d.model+"<br />Start-End: "+start.toString()+"-"+end.toString()+"cM<br />Range: "+range.toString()+"cM");
}


var gen_scatter_size = function(d) {
    return(d.size);
}

var drawCircos = function (error, karyotypes, layout, trait_files, lod_files, scatter_configs, stack_configs, line_configs, circos_trait_config) {
  if(error) {
    throw error;
  }
  circosScatter
    .layout(
        karyotypes,
        layout
    );
  for( i = 0; i < circos_trait_config.length; i++ ) {
    circos_trait2traitname[circos_trait_config[i].trait] = circos_trait_config[i].label;
  }
  var q = d3.queue();
  for( i = 0; i < trait_files.length; i++ ) {
    trait_file                =  gfolder_prefix + trait_files[i];
    scatter_configs[i].color          =  gen_scatter_color;
    scatter_configs[i].tooltipContent =  gen_scatter_tooltip;
    scatter_configs[i].size           =  gen_scatter_size;
    for( j = 0; j < stack_configs[i].length; j++ ) {
        stack_configs[i][j].color          =  gen_scatter_color;
        stack_configs[i][j].tooltipContent =  gen_stack_tooltip;
    }
    q.defer(d3.csv, trait_file)
     .defer(get_current_config, scatter_configs[i])
     .defer(get_current_config, stack_configs[i])
  }
  q.awaitAll(inject_scatter);

  var q = d3.queue();
  for( i = 0; i < lod_files.length; i++ ) {
    lod_file = gfolder_prefix + lod_files[i]
    q.defer(d3.json, lod_file)
     .defer(get_current_config, line_configs[i])
  }
  q.awaitAll(inject_line);
}



var generate_new_circos = function (container="#scatterChart", folder_prefix="", width="1024", height="1024", opacity=0.75, highlight_opacity=1.0, bubble_scale=15) {
    circosScatter = new Circos({
        container: container,
        width: width,
        height: height 
    });

    gfolder_prefix                = folder_prefix
    qtl_bubble_opacity            = opacity;
    qtl_bubble_highlight_opacity  = highlight_opacity;
    qtl_bubble_scale_factor       = bubble_scale;

    d3.queue()
    .defer(d3.json, folder_prefix+'karyotype.json')
    .defer(d3.json, folder_prefix+'configs/circos/layout.json')
    .defer(d3.json, folder_prefix+'configs/circos/all_traits.json')
    .defer(d3.json, folder_prefix+'configs/circos/all_lods.json')
    .defer(d3.json, folder_prefix+'configs/circos/scatter.configs.json')
    .defer(d3.json, folder_prefix+'configs/circos/stack.configs.json')
    .defer(d3.json, folder_prefix+'configs/circos/line.configs.json')
    .defer(d3.csv,  folder_prefix+'configs/model-traits.cfg.csv')
    .await(drawCircos);
}
