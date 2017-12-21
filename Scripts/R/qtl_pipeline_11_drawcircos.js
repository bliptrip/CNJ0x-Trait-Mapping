var width = document.getElementsByClassName('mdl-card__supporting-text')[0].offsetWidth
var circosScatter = new Circos({
      container: '#scatterChart',
      width: width,
      height: width
});

var get_current_config = function(config,callback) {
    setTimeout(function() {
        callback(null, config);
    }, 250);
}

var show_lod_profs = function(datum, index, nodes, event) {
    var mclass = 'lines--'+datum.model+'--'+datum.mtraits+'--'+datum.trait;
    //Hide other lodprofiles that are displaying in the current trait
    d3.selectAll("*[class^='lines--'][class$='"+datum.mtraits+"--"+datum.trait+"']").attr('visibility','hidden');
    d3.selectAll("*[class='"+mclass+"']").attr('visibility','visible');
}

var hide_lod_profs = function(datum, index, nodes, event) {
    var mclass = 'lines--'+datum[0].model+'--'+datum[0].mtraits+'--'+datum[0].trait;
    d3.selectAll("*[class='"+mclass+"']").attr('visibility','hidden');
}

var highlight_markers = function(datum, index, nodes, event) {
    var marker = datum.nearest_marker;
    d3.selectAll("*[marker='"+marker+"']").attr('opacity', '1.0');
}

var unhighlight_markers = function(datum, index, nodes, event) {
    var marker = datum.nearest_marker;
    d3.selectAll("*[marker='"+marker+"']").attr('opacity', '0.75');
}

var inject_scatter = function(error, data)
{
    bin_size = 3; //Specifies the number of indices per 'bin' of data
    for( i = 0; i < (data.length/bin_size); i++ ) {
        trait_data = data[bin_size*i];
        scatter_config     = data[(bin_size*i)+1];
        scatter_config.events = {
            'click.showlodprofs': show_lod_profs

        };
        scatter_trait_data = trait_data.map(function(d) {
                return({
                    block_id: "vm"+d.chr,
                    position: +d.position*1000,
                    value: +d.model_idx,  
                    model: d.model,
                    mtraits: d.mtraits,
                    trait: d.trait,
                    color: d.color,
                    stroke_color: d.stroke_color,
                    model_color: d.model_color,
                    nearest_marker: d.nearest_marker,
                    position: d.position * 1000,
                    size: Math.round(d.marker_variance * 10.0),
                    class: d.class
                });
            });
        mtraits_trait = "scatter--" + scatter_trait_data[0].mtraits + "--" + scatter_trait_data[0].trait;
        circosScatter.scatter(mtraits_trait, scatter_trait_data, scatter_config);
        stack_configs     = data[(bin_size*i)+2];
        for( j = 0; j < stack_configs.length; j++ ) {
            stack_trait_data = trait_data.filter(function(d) 
                    { 
                        var match = (+d.model_idx == (j+1)); 
                        return(match);
                    })
                    .map(function(d) {
                        return({
                            block_id: "vm"+d.chr,
                            start: (+d.position - (+d.interval/2))*1000,
                            end: (+d.position + (+d.interval/2))*1000,
                            model: d.model,
                            mtraits: d.mtraits,
                            trait: d.trait,
                            color: d.color,
                            stroke_color: d.stroke_color,
                            model_color: d.model_color,
                            nearest_marker: d.nearest_marker,
                            class: d.class
                        });
                    });
            if( stack_trait_data.length > 0 ) {
                mtraits_trait = "stacks--"+stack_trait_data[0].model+"--"+ stack_trait_data[0].mtraits+"--"+stack_trait_data[0].trait;
                circosScatter.stack(mtraits_trait, stack_trait_data, stack_configs[j]);
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
            lod_trait_data = lod_data[model].map(function(d) {
            //I need to figure out how to separate the different models into different datasets
                    return({
                        block_id: "vm"+d.chr,
                        position: +d.position*1000,
                        value: +d.lod,
                        model: d.model,
                        mtraits: d.mtraits,
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
            mtraits_trait = "lines--"+model+"--"+lod_trait_data[0].mtraits+"--"+lod_trait_data[0].trait;
            circosScatter.line(mtraits_trait, lod_trait_data, lod_config);
        }
    }
    circosScatter.render();
    d3.selectAll(".point")
      .on('mouseover.highlight', highlight_markers)
      .on('mouseout.highlight', unhighlight_markers);
    d3.selectAll(".tile")
      .on('mouseover.highlight', highlight_markers)
      .on('mouseout.highlight', unhighlight_markers);
    d3.selectAll("*[class^='lines--']").attr('visibility','hidden');
    d3.selectAll("g[class='line']")
      .on('click', hide_lod_profs);
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
    return(d.model+"<br>"+d.trait+"<br>"+d.nearest_marker);
}

var gen_scatter_size = function(d) {
    return(d.size);
}

var drawCircos = function (error, karyotypes, layout, trait_files, lod_files, scatter_configs, stack_configs, line_configs) {
  if(error) {
    throw error;
  }
  circosScatter
    .layout(
        karyotypes,
        layout
    );
  var q = d3.queue();
  for( i = 0; i < trait_files.length; i++ ) {
    trait_file                =  trait_files[i];
    scatter_configs[i].color          =  gen_scatter_color;
    //scatter_configs[i].strokeColor    =  gen_scatter_stroke_color;
    scatter_configs[i].tooltipContent =  gen_scatter_tooltip;
    scatter_configs[i].size           =  gen_scatter_size;
    for( j = 0; j < stack_configs[i].length; j++ ) {
        stack_configs[i][j].color          =  gen_scatter_color;
        //scatter_configs[i].strokeColor    =  gen_scatter_stroke_color;
        stack_configs[i][j].tooltipContent =  gen_scatter_tooltip;
    }
    q.defer(d3.csv, trait_file)
     .defer(get_current_config, scatter_configs[i])
     .defer(get_current_config, stack_configs[i])
  }
  q.awaitAll(inject_scatter);

  var q = d3.queue();
  for( i = 0; i < lod_files.length; i++ ) {
    q.defer(d3.json, lod_files[i])
     .defer(get_current_config, line_configs[i])
  }
  q.awaitAll(inject_line);
}

d3.queue()
  .defer(d3.json, '../../Workflows/1/traits/plots/circos/karyotype.json')
  .defer(d3.json, '../../Workflows/1/configs/circos/layout.json')
  .defer(d3.json, '../../Workflows/1/configs/circos/all_traits.json')
  .defer(d3.json, '../../Workflows/1/configs/circos/all_lods.json')
  .defer(d3.json, '../../Workflows/1/configs/circos/scatter.configs.json')
  .defer(d3.json, '../../Workflows/1/configs/circos/stack.configs.json')
  .defer(d3.json, '../../Workflows/1/configs/circos/line.configs.json')
  .await(drawCircos);
