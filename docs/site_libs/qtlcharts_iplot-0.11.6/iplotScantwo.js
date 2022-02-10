"use strict";

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

// Generated by CoffeeScript 2.2.2
// iplotScantwo: interactive plot of scantwo results (2-dim, 2-QTL genome scan)
// Karl W Broman
var add_symmetric_lod, iplotScantwo, lod_for_heatmap;

iplotScantwo = function iplotScantwo(widgetdiv, scantwo_data, pheno_and_geno, chartOpts) {
  var add_cell_tooltips, altrectcolor, axispos, boxcolor, boxwidth, chrGap, chrlinecolor, chrlinewidth, cicolors, color, cur_chr1, cur_chr2, div, eff_hpos, eff_vpos, form, g_eff, g_heatmap, g_scans, gn, hbot, heatmap_height, heatmap_width, height, hright, i, left, leftsel, leftvalue, linecolor, linewidth, margin, mycichart, mydotchart, mylod2dheatmap, mylodchart, n, ncat, nullcolor, nyticks_ci, nyticks_lod, nyticks_pxg, oneAtTop, options, plot_effects, plot_scan, pointsize, pointstroke, rectcolor, ref, ref1, ref10, ref11, ref12, ref13, ref14, ref15, ref16, ref17, ref18, ref19, ref2, ref20, ref21, ref22, ref23, ref24, ref25, ref26, ref27, ref28, ref29, ref3, ref30, ref31, ref32, ref33, ref4, ref5, ref6, ref7, ref8, ref9, right, rightsel, rightvalue, scans_hpos, scans_vpos, segstrokewidth, segwidth, submit, svg, titlepos, w, wbot, widgetdivid, width, wright, x, xlab_lod, ylab_eff, ylab_lod, yticks_ci, yticks_lod, yticks_pxg, zthresh;
  // chartOpts start
  height = (ref = chartOpts != null ? chartOpts.height : void 0) != null ? ref : 1200; // total height of chart in pixels
  width = (ref1 = chartOpts != null ? chartOpts.width : void 0) != null ? ref1 : 1100; // total width of chart in pixels
  chrGap = (ref2 = chartOpts != null ? chartOpts.chrGap : void 0) != null ? ref2 : 2; // gaps between chr in heat map
  wright = (ref3 = chartOpts != null ? chartOpts.wright : void 0) != null ? ref3 : width / 2; // width (in pixels) of right panels
  hbot = (ref4 = chartOpts != null ? chartOpts.hbot : void 0) != null ? ref4 : height / 5; // height (in pixels) of each of the lower panels
  margin = (ref5 = chartOpts != null ? chartOpts.margin : void 0) != null ? ref5 : {
    left: 60,
    top: 50,
    right: 10,
    bottom: 40,
    inner: 5 // margins in each panel
  };
  axispos = (ref6 = chartOpts != null ? chartOpts.axispos : void 0) != null ? ref6 : {
    xtitle: 25,
    ytitle: 30,
    xlabel: 5,
    ylabel: 5 // axis positions in heatmap
  };
  titlepos = (ref7 = chartOpts != null ? chartOpts.titlepos : void 0) != null ? ref7 : 20; // position of chart title in pixels
  rectcolor = (ref8 = chartOpts != null ? chartOpts.rectcolor : void 0) != null ? ref8 : "#e6e6e6"; // color for background rectangle
  altrectcolor = (ref9 = chartOpts != null ? chartOpts.altrectcolor : void 0) != null ? ref9 : "#c8c8c8"; // alternate rectangle in lower panels
  chrlinecolor = (ref10 = chartOpts != null ? chartOpts.chrlinecolor : void 0) != null ? ref10 : ""; // color of lines between chromosomes (if "", leave off)
  chrlinewidth = (ref11 = chartOpts != null ? chartOpts.chrlinewidth : void 0) != null ? ref11 : 2; // width of lines between chromosomes
  nullcolor = (ref12 = chartOpts != null ? chartOpts.nullcolor : void 0) != null ? ref12 : "#e6e6e6"; // color of null pixels in heat map
  boxcolor = (ref13 = chartOpts != null ? chartOpts.boxcolor : void 0) != null ? ref13 : "black"; // color of box around each panel
  boxwidth = (ref14 = chartOpts != null ? chartOpts.boxwidth : void 0) != null ? ref14 : 2; // width of box around each panel
  linecolor = (ref15 = chartOpts != null ? chartOpts.linecolor : void 0) != null ? ref15 : "slateblue"; // line color in lower panels
  linewidth = (ref16 = chartOpts != null ? chartOpts.linewidth : void 0) != null ? ref16 : 2; // line width in lower panels
  pointsize = (ref17 = chartOpts != null ? chartOpts.pointsize : void 0) != null ? ref17 : 2; // point size in right panels
  pointstroke = (ref18 = chartOpts != null ? chartOpts.pointstroke : void 0) != null ? ref18 : "black"; // color of outer circle in right panels
  cicolors = (ref19 = chartOpts != null ? chartOpts.cicolors : void 0) != null ? ref19 : null; // colors for CIs in QTL effect plot; also used for points in phe x gen plot
  segwidth = (ref20 = chartOpts != null ? chartOpts.segwidth : void 0) != null ? ref20 : 0.4; // segment width in CI chart as proportion of distance between categories
  segstrokewidth = (ref21 = chartOpts != null ? chartOpts.segstrokewidth : void 0) != null ? ref21 : 3; // stroke width for segments in CI chart
  color = (ref22 = chartOpts != null ? chartOpts.color : void 0) != null ? ref22 : "slateblue"; // color for heat map
  oneAtTop = (ref23 = chartOpts != null ? chartOpts.oneAtTop : void 0) != null ? ref23 : false; // whether to put chr 1 at top of heatmap
  zthresh = (ref24 = chartOpts != null ? chartOpts.zthresh : void 0) != null ? ref24 : 0; // LOD values below this threshold aren't shown (on LOD_full scale)
  ylab_eff = (ref25 = chartOpts != null ? chartOpts.ylab_eff : void 0) != null ? ref25 : "Phenotype"; // y-axis label in dot and ci charts
  xlab_lod = (ref26 = chartOpts != null ? chartOpts.xlab_lod : void 0) != null ? ref26 : "Chromosome"; // x-axis label in lod charts
  ylab_lod = (ref27 = chartOpts != null ? chartOpts.ylab_lod : void 0) != null ? ref27 : "LOD score"; // y-axis label in lod charts
  nyticks_lod = (ref28 = chartOpts != null ? chartOpts.nyticks_lod : void 0) != null ? ref28 : 5; // no. ticks on y-axis in LOD curve panels
  yticks_lod = (ref29 = chartOpts != null ? chartOpts.yticks_lod : void 0) != null ? ref29 : null; // vector of tick positions on y-axis in LOD curve panels
  nyticks_ci = (ref30 = chartOpts != null ? chartOpts.nyticks_ci : void 0) != null ? ref30 : 5; // no. ticks on y-axis in CI panel
  yticks_ci = (ref31 = chartOpts != null ? chartOpts.yticks_ci : void 0) != null ? ref31 : null; // vector of tick positions on y-axis in CI panel
  nyticks_pxg = (ref32 = chartOpts != null ? chartOpts.nyticks_pxg : void 0) != null ? ref32 : 5; // no. ticks on y-axis in dot chart of phenotype x genotype
  yticks_pxg = (ref33 = chartOpts != null ? chartOpts.yticks_pxg : void 0) != null ? ref33 : null; // vector of tick positions on y-axis in dot chart of phenotype x genotype
  // chartOpts end

  // make sure list args have all necessary bits
  margin = d3panels.check_listarg_v_default(margin, {
    left: 60,
    top: 50,
    right: 10,
    bottom: 40,
    inner: 5
  });
  axispos = d3panels.check_listarg_v_default(axispos, {
    xtitle: 25,
    ytitle: 30,
    xlabel: 5,
    ylabel: 5
  });
  // htmlwidget div element containing the chart, and its ID
  div = d3.select(widgetdiv);
  widgetdivid = div.attr("id");
  svg = div.select("svg");
  // force chrnames to be a list
  scantwo_data.chrnames = d3panels.forceAsArray(scantwo_data.chrnames);
  scantwo_data.nmar = d3panels.forceAsArray(scantwo_data.nmar);
  // size of heatmap region
  w = d3.min([height - hbot * 2, width - wright]);
  heatmap_width = w;
  heatmap_height = w;
  hright = heatmap_height / 2;
  width = heatmap_width + wright;
  height = heatmap_height + hbot * 2;
  wbot = width / 2;
  // selected LODs on left and right
  leftvalue = "int";
  rightvalue = "fv1";
  // keep track of chromosome heatmap selections
  cur_chr1 = cur_chr2 = '';
  // cicolors: check they're the right length or construct them
  if (pheno_and_geno != null) {
    gn = pheno_and_geno.genonames;
    ncat = d3.max(function () {
      var results;
      results = [];
      for (x in gn) {
        results.push(gn[x].length);
      }
      return results;
    }());
    if (cicolors != null) {
      cicolors = d3panels.expand2vector(cicolors, ncat);
      n = cicolors.length;
      if (n < ncat) {
        // not enough, display error
        d3panels.displayError("length(cicolors) (" + n + ") < maximum no. genotypes (" + ncat + ")");
        cicolors = function () {
          var l, ref34, results;
          // not provided; select them
          results = [];
          for (i = l = 0, ref34 = ncat; 0 <= ref34 ? l < ref34 : l > ref34; i = 0 <= ref34 ? ++l : --l) {
            results.push(cicolors[i % n]);
          }
          return results;
        }();
      }
    } else {
      cicolors = d3panels.selectGroupColors(ncat, "dark"); // cicolors provided; expand to ncat
    }
  }
  // drop-down menus
  options = ["full", "fv1", "int", "add", "av1"];
  form = div.insert("div", ":first-child").attr("id", "form").attr("class", "qtlcharts").attr("height", "24px");
  left = form.append("div").text(oneAtTop ? "bottom-left: " : "top-left: ").style("float", "left").style("margin-left", "50px");
  leftsel = left.append("select").attr("id", "leftselect_" + widgetdivid).attr("name", "left");
  leftsel.selectAll("empty").data(options).enter().append("option").attr("value", function (d) {
    return d;
  }).text(function (d) {
    return d;
  }).attr("selected", function (d) {
    if (d === leftvalue) {
      return "selected";
    }
    return null;
  });
  right = form.append("div").text(oneAtTop ? "top-right: " : "bottom-right: ").style("float", "left").style("margin-left", "50px");
  rightsel = right.append("select").attr("id", "rightselect_" + widgetdivid).attr("name", "right");
  rightsel.selectAll("empty").data(options).enter().append("option").attr("value", function (d) {
    return d;
  }).text(function (d) {
    return d;
  }).attr("selected", function (d) {
    if (d === rightvalue) {
      return "selected";
    }
    return null;
  });
  submit = form.append("div").style("float", "left").style("margin-left", "50px").append("button").attr("name", "refresh").text("Refresh").on("click", function () {
    cur_chr1 = cur_chr2 = '';
    leftsel = document.getElementById("leftselect_" + widgetdivid);
    leftvalue = leftsel.options[leftsel.selectedIndex].value;
    rightsel = document.getElementById("rightselect_" + widgetdivid);
    rightvalue = rightsel.options[rightsel.selectedIndex].value;
    scantwo_data.lod = lod_for_heatmap(scantwo_data, leftvalue, rightvalue);
    mylod2dheatmap.remove();
    mylod2dheatmap(div.select("g#chrheatmap"), scantwo_data);
    return add_cell_tooltips();
  });
  // add the full,add,int,fv1,av1 lod matrices to scantwo_data
  // (and remove the non-symmetric ones)
  scantwo_data = add_symmetric_lod(scantwo_data);
  scantwo_data.lod = lod_for_heatmap(scantwo_data, leftvalue, rightvalue);
  mylod2dheatmap = d3panels.lod2dheatmap({
    height: heatmap_height,
    width: heatmap_width,
    margin: margin,
    axispos: axispos,
    chrGap: chrGap,
    chrlinecolor: chrlinecolor,
    chrlinewidth: chrlinewidth,
    xlab: xlab_lod,
    ylab: ylab_lod,
    rectcolor: "white",
    nullcolor: nullcolor,
    boxcolor: boxcolor,
    boxwidth: boxwidth,
    colors: ["white", color],
    zlim: [0, scantwo_data.max.full],
    zthresh: zthresh,
    oneAtTop: oneAtTop,
    tipclass: widgetdivid
  });
  g_heatmap = svg.append("g").attr("id", "chrheatmap");
  mylod2dheatmap(g_heatmap, scantwo_data);
  // function to add tool tips and handle clicking
  add_cell_tooltips = function add_cell_tooltips() {
    mylod2dheatmap.celltip().html(function (d) {
      var leftlod, mari, marj, rightlod;
      mari = scantwo_data.marker[d.xindex];
      marj = scantwo_data.marker[d.yindex];
      if (+d.xindex > +d.yindex) {
        // +'s ensure number not string
        leftlod = d3.format(".1f")(scantwo_data[leftvalue][d.xindex][d.yindex]);
        rightlod = d3.format(".1f")(scantwo_data[rightvalue][d.yindex][d.xindex]);
        return "(" + marj + " " + mari + ") " + rightvalue + " = " + rightlod + ", " + leftvalue + " = " + leftlod;
      } else if (+d.yindex > +d.xindex) {
        leftlod = d3.format(".1f")(scantwo_data[leftvalue][d.yindex][d.xindex]);
        rightlod = d3.format(".1f")(scantwo_data[rightvalue][d.xindex][d.yindex]);
        return "(" + marj + " " + mari + ") " + leftvalue + " = " + leftlod + ", " + rightvalue + " = " + rightlod;
      } else {
        return mari;
      }
    });
    return mylod2dheatmap.cells().on("click", function (d) {
      var mari, marj;
      mari = scantwo_data.marker[d.xindex];
      marj = scantwo_data.marker[d.yindex];
      if (d.xindex === d.yindex) {
        // skip the diagonal case
        return null;
      }
      // plot the cross-sections as genome scans, below
      plot_scan(d.xindex, 0, 0, leftvalue);
      plot_scan(d.xindex, 1, 0, rightvalue);
      plot_scan(d.yindex, 0, 1, leftvalue);
      plot_scan(d.yindex, 1, 1, rightvalue);
      // plot the effect plot and phe x gen plot to right
      if (pheno_and_geno != null) {
        return plot_effects(d.xindex, d.yindex);
      }
    });
  };
  add_cell_tooltips();
  // to hold groups and positions of scan and effect plots
  mylodchart = [[null, null], [null, null]];
  scans_hpos = [0, wbot];
  scans_vpos = [heatmap_height, heatmap_height + hbot];
  mydotchart = null;
  mycichart = null;
  eff_hpos = [heatmap_width, heatmap_width];
  eff_vpos = [0, heatmap_height / 2];
  g_scans = [[null, null], [null, null]];
  plot_scan = function plot_scan(markerindex, panelrow, panelcol, lod) {
    var data;
    data = {
      chrname: scantwo_data.chrnames,
      chr: scantwo_data.chr,
      pos: scantwo_data.pos,
      lod: function () {
        var l, len, ref34, results;
        ref34 = scantwo_data[lod][markerindex];
        results = [];
        for (l = 0, len = ref34.length; l < len; l++) {
          x = ref34[l];
          results.push(x);
        }
        return results;
      }(),
      marker: scantwo_data.marker
    };
    if (mylodchart[panelrow][panelcol] != null) {
      mylodchart[panelrow][panelcol].remove();
    }
    mylodchart[panelrow][panelcol] = d3panels.lodchart({
      height: hbot,
      width: wbot,
      margin: margin,
      axispos: axispos,
      ylim: [0.0, scantwo_data.max[lod] * 1.05],
      nyticks: nyticks_lod,
      yticks: yticks_lod,
      rectcolor: rectcolor,
      altrectcolor: altrectcolor,
      chrlinecolor: chrlinecolor,
      chrlinewidth: chrlinewidth,
      boxcolor: boxcolor,
      boxwidth: boxwidth,
      linewidth: linewidth,
      linecolor: linecolor,
      pointsize: 0,
      pointcolor: "",
      pointstroke: "",
      lodvarname: "lod",
      chrGap: chrGap,
      xlab: xlab_lod,
      ylab: ylab_lod,
      title: data.marker[markerindex] + " : " + lod,
      titlepos: titlepos,
      tipclass: widgetdivid
    });
    if (g_scans[panelrow][panelcol] == null) {
      g_scans[panelrow][panelcol] = svg.append("g").attr("id", "scan_" + (panelrow + 1) + "_" + (panelcol + 1)).attr("transform", "translate(" + scans_hpos[panelcol] + ", " + scans_vpos[panelrow] + ")");
    }
    return mylodchart[panelrow][panelcol](g_scans[panelrow][panelcol], data);
  };
  g_eff = [null, null];
  plot_effects = function plot_effects(markerindex1, markerindex2) {
    var _d3panels$cichart;

    var chr1, chr2, chrtype1, chrtype2, ci_data, cicolors_expanded, cis, dpos, effcharts, fgnames, force, g, g1, g2, gn1, gn2, gnames1, gnames2, j, k, l, m, mar1, mar2, mgnames, ng1, ng2, ngf, ngm, o, p, point_jitter, points, pos1, pxg_data, q, r, ref34, ref35, ref36, ref37, ref38, ref39, ref40, ref41, ref42, results, s, scaledPoints, t, tmp, v, xscale;
    mar1 = scantwo_data.marker[markerindex1];
    mar2 = scantwo_data.marker[markerindex2];
    g1 = pheno_and_geno.geno[mar1];
    g2 = pheno_and_geno.geno[mar2];
    chr1 = pheno_and_geno.chr[mar1];
    chr2 = pheno_and_geno.chr[mar2];
    chrtype1 = pheno_and_geno.chrtype[chr1];
    chrtype2 = pheno_and_geno.chrtype[chr2];
    g = [];
    gn1 = [];
    gn2 = [];
    cicolors_expanded = [];
    // need to deal separately with X chr
    // [this mess is because if females are AA/AB/BB and males AY/BY
    //  we want to just show 3x3 + 2x2 = 13 possible two-locus genotypes,
    //  not all (3+2)x(3+2) = 25]
    if (chr1 === chr2 && chrtype1 === "X" && pheno_and_geno.X_geno_by_sex != null) {
      fgnames = pheno_and_geno.X_geno_by_sex[0];
      mgnames = pheno_and_geno.X_geno_by_sex[1];
      ngf = fgnames.length;
      ngm = mgnames.length;
      tmp = function () {
        var results = [];
        for (var l = 0, ref34 = ngf + ngm; 0 <= ref34 ? l < ref34 : l > ref34; 0 <= ref34 ? l++ : l--) {
          results.push(l);
        }
        return results;
      }.apply(this);
      m = function () {
        var results;
        results = [];
        for (j in tmp) {
          results.push(function () {
            var results1;
            results1 = [];
            for (i in tmp) {
              results1.push(-1);
            }
            return results1;
          }());
        }
        return results;
      }();
      k = 0;
      for (i = l = 0, ref35 = ngf; 0 <= ref35 ? l < ref35 : l > ref35; i = 0 <= ref35 ? ++l : --l) {
        for (j = o = 0, ref36 = ngf; 0 <= ref36 ? o < ref36 : o > ref36; j = 0 <= ref36 ? ++o : --o) {
          gn1.push(fgnames[j]);
          gn2.push(fgnames[i]);
          cicolors_expanded.push(cicolors[i]);
          m[i][j] = k;
          k++;
        }
      }
      for (i = q = 0, ref37 = ngm; 0 <= ref37 ? q < ref37 : q > ref37; i = 0 <= ref37 ? ++q : --q) {
        for (j = r = 0, ref38 = ngm; 0 <= ref38 ? r < ref38 : r > ref38; j = 0 <= ref38 ? ++r : --r) {
          gn1.push(mgnames[j]);
          gn2.push(mgnames[i]);
          cicolors_expanded.push(cicolors[i]);
          m[i + ngf][j + ngf] = k;
          k++;
        }
      }
      g = function () {
        var results;
        results = [];
        for (i in g1) {
          results.push(m[g1[i] - 1][g2[i] - 1] + 1);
        }
        return results;
      }();
    } else {
      gnames1 = pheno_and_geno.genonames[chr1];
      gnames2 = pheno_and_geno.genonames[chr2];
      ng1 = gnames1.length;
      ng2 = gnames2.length;
      g = function () {
        var results;
        results = [];
        for (i in g1) {
          results.push(g1[i] + (g2[i] - 1) * ng1);
        }
        return results;
      }();
      for (i = s = 0, ref39 = ng2; 0 <= ref39 ? s < ref39 : s > ref39; i = 0 <= ref39 ? ++s : --s) {
        for (j = t = 0, ref40 = ng1; 0 <= ref40 ? t < ref40 : t > ref40; j = 0 <= ref40 ? ++t : --t) {
          gn1.push(gnames1[j]);
          gn2.push(gnames2[i]);
          cicolors_expanded.push(cicolors[i]);
        }
      }
    }
    pxg_data = {
      x: g,
      y: pheno_and_geno.pheno,
      indID: pheno_and_geno.indID
    };
    if (mycichart != null) {
      // remove the CI chart no matter what
      mycichart.remove();
    }
    if (cur_chr1 !== chr1 || cur_chr2 !== chr2) {
      if (mydotchart != null) {
        mydotchart.remove();
      }
      mydotchart = d3panels.dotchart({
        height: hright,
        width: wright,
        margin: margin,
        axispos: axispos,
        rectcolor: rectcolor,
        boxcolor: boxcolor,
        boxwidth: boxwidth,
        pointsize: pointsize,
        pointstroke: pointstroke,
        xcategories: function () {
          var results = [];
          for (var v = 1, ref41 = gn1.length; 1 <= ref41 ? v <= ref41 : v >= ref41; 1 <= ref41 ? v++ : v--) {
            results.push(v);
          }
          return results;
        }.apply(this),
        xcatlabels: gn1,
        xlab: "",
        ylab: ylab_eff,
        nyticks: nyticks_pxg,
        yticks: yticks_pxg,
        dataByInd: false,
        title: mar1 + " : " + mar2,
        titlepos: titlepos,
        tipclass: widgetdivid
      });
      if (g_eff[1] == null) {
        g_eff[1] = svg.append("g").attr("id", "eff_1").attr("transform", "translate(" + eff_hpos[1] + ", " + eff_vpos[1] + ")");
      }
      mydotchart(g_eff[1], pxg_data);
      // revise point colors
      mydotchart.points().attr("fill", function (d, i) {
        return cicolors_expanded[g[i] - 1];
      });
    } else {
      // remove marker text
      d3.select("#markerlab1").remove();
      d3.select("#xaxislab1").remove();
      // grab scale and get info to take inverse
      xscale = mydotchart.xscale();
      pos1 = xscale(1);
      dpos = xscale(2) - xscale(1);
      point_jitter = function point_jitter(d) {
        var u;
        u = (d - pos1) / dpos;
        return u - Math.round(u);
      };
      // move points to new x-axis position
      points = mydotchart.points().transition().duration(1000).attr("cx", function (d, i) {
        var cx, u;
        cx = d3.select(this).attr("cx");
        u = point_jitter(cx);
        return xscale(g[i] + u);
      }).attr("fill", function (d, i) {
        return cicolors_expanded[g[i] - 1];
      });
      // use force to move them apart again
      scaledPoints = [];
      points.each(function (d, i) {
        return scaledPoints.push({
          x: +d3.select(this).attr("cx"),
          y: +d3.select(this).attr("cy"),
          fy: +d3.select(this).attr("cy"),
          truex: xscale(g[i])
        });
      });
      force = d3.forceSimulation(scaledPoints).force("x", d3.forceX(function (d) {
        return d.truex;
      })).force("collide", d3.forceCollide(pointsize * 1.1)).stop();
      (function () {
        var results = [];
        for (var v = 0; v <= 30; v++) {
          results.push(v);
        }
        return results;
      }).apply(this).map(function (d) {
        force.tick();
        return points.attr("cx", function (d, i) {
          return scaledPoints[i].x;
        });
      });
    }
    cur_chr1 = chr1;
    cur_chr2 = chr2;
    cis = d3panels.ci_by_group(g, pheno_and_geno.pheno, 2);
    ci_data = {
      mean: function () {
        var ref42, ref43, ref44, results, v;
        results = [];
        for (x = v = 1, ref42 = gn1.length; 1 <= ref42 ? v <= ref42 : v >= ref42; x = 1 <= ref42 ? ++v : --v) {
          results.push((ref43 = (ref44 = cis[x]) != null ? ref44.mean : void 0) != null ? ref43 : null);
        }
        return results;
      }(),
      low: function () {
        var ref42, ref43, ref44, results, v;
        results = [];
        for (x = v = 1, ref42 = gn1.length; 1 <= ref42 ? v <= ref42 : v >= ref42; x = 1 <= ref42 ? ++v : --v) {
          results.push((ref43 = (ref44 = cis[x]) != null ? ref44.low : void 0) != null ? ref43 : null);
        }
        return results;
      }(),
      high: function () {
        var ref42, ref43, ref44, results, v;
        results = [];
        for (x = v = 1, ref42 = gn1.length; 1 <= ref42 ? v <= ref42 : v >= ref42; x = 1 <= ref42 ? ++v : --v) {
          results.push((ref43 = (ref44 = cis[x]) != null ? ref44.high : void 0) != null ? ref43 : null);
        }
        return results;
      }(),
      categories: function () {
        var results = [];
        for (var v = 1, ref42 = gn1.length; 1 <= ref42 ? v <= ref42 : v >= ref42; 1 <= ref42 ? v++ : v--) {
          results.push(v);
        }
        return results;
      }.apply(this)
    };
    mycichart = d3panels.cichart((_d3panels$cichart = {
      height: hright,
      width: wright,
      margin: margin,
      axispos: axispos,
      rectcolor: rectcolor,
      boxcolor: boxcolor,
      boxwidth: boxwidth,
      segcolor: cicolors_expanded,
      segwidth: segwidth,
      segstrokewidth: segstrokewidth,
      vertsegcolor: cicolors_expanded
    }, _defineProperty(_d3panels$cichart, "segstrokewidth", linewidth), _defineProperty(_d3panels$cichart, "xlab", ""), _defineProperty(_d3panels$cichart, "ylab", ylab_eff), _defineProperty(_d3panels$cichart, "nyticks", nyticks_ci), _defineProperty(_d3panels$cichart, "yticks", yticks_ci), _defineProperty(_d3panels$cichart, "xcatlabels", gn1), _defineProperty(_d3panels$cichart, "title", mar1 + " : " + mar2), _defineProperty(_d3panels$cichart, "titlepos", titlepos), _defineProperty(_d3panels$cichart, "tipclass", widgetdivid), _d3panels$cichart));
    if (g_eff[0] == null) {
      g_eff[0] = svg.append("g").attr("id", "eff_0").attr("transform", "translate(" + eff_hpos[0] + ", " + eff_vpos[0] + ")");
    }
    mycichart(g_eff[0], ci_data);
    effcharts = [mycichart, mydotchart];
    // add second row of labels
    results = [];
    for (p = v = 0; v <= 1; p = ++v) {
      effcharts[p].svg().append("g").attr("class", "x axis").attr("id", "xaxislab" + p // second row of genotypes
      ).selectAll("empty").data(gn2).enter().append("text").attr("x", function (d, i) {
        return mydotchart.xscale()(i + 1);
      }).attr("y", hright - margin.bottom / 2 + axispos.xlabel).text(function (d) {
        return d;
      });
      results.push(effcharts[p].svg().append("g").attr("class", "x axis").attr("id", "markerlab" + p // marker name labels
      ).selectAll("empty").data([mar1, mar2]).enter().append("text").attr("x", (margin.left + mydotchart.xscale()(1)) / 2.0).attr("y", function (d, i) {
        return hright - margin.bottom / (i + 1) + axispos.xlabel;
      }).style("text-anchor", "end").text(function (d) {
        return d + ":";
      }));
    }
    return results;
  };
  if (chartOpts.heading != null) {
    d3.select("div#htmlwidget_container").insert("h2", ":first-child").html(chartOpts.heading).style("font-family", "sans-serif");
  }
  if (chartOpts.caption != null) {
    d3.select("body").append("p").attr("class", "caption").html(chartOpts.caption);
  }
  if (chartOpts.footer != null) {
    return d3.select("body").append("div").html(chartOpts.footer).style("font-family", "sans-serif");
  }
};

// add full,add,int,av1,fv1 lod scores to scantwo_data
add_symmetric_lod = function add_symmetric_lod(scantwo_data) {
  var i, j, l, len, o, q, r, ref, ref1, ref2, ref3, ref4, ref5, s;
  scantwo_data.full = scantwo_data.lod.map(function (d) {
    return d.map(function (dd) {
      return dd;
    });
  });
  scantwo_data.add = scantwo_data.lod.map(function (d) {
    return d.map(function (dd) {
      return dd;
    });
  });
  scantwo_data.fv1 = scantwo_data.lodv1.map(function (d) {
    return d.map(function (dd) {
      return dd;
    });
  });
  scantwo_data.av1 = scantwo_data.lodv1.map(function (d) {
    return d.map(function (dd) {
      return dd;
    });
  });
  scantwo_data.int = scantwo_data.lod.map(function (d) {
    return d.map(function (dd) {
      return dd;
    });
  });
  for (i = l = 0, ref = scantwo_data.lod.length - 1; 0 <= ref ? l < ref : l > ref; i = 0 <= ref ? ++l : --l) {
    for (j = o = ref1 = i, ref2 = scantwo_data.lod[i].length; ref1 <= ref2 ? o < ref2 : o > ref2; j = ref1 <= ref2 ? ++o : --o) {
      scantwo_data.full[i][j] = scantwo_data.lod[j][i];
      scantwo_data.add[j][i] = scantwo_data.lod[i][j];
      scantwo_data.fv1[i][j] = scantwo_data.lodv1[j][i];
      scantwo_data.av1[j][i] = scantwo_data.lodv1[i][j];
    }
  }
  scantwo_data.one = [];
  for (i = q = 0, ref3 = scantwo_data.lod.length; 0 <= ref3 ? q < ref3 : q > ref3; i = 0 <= ref3 ? ++q : --q) {
    scantwo_data.one.push(scantwo_data.lod[i]);
    for (j = r = 0, ref4 = scantwo_data.lod.length; 0 <= ref4 ? r < ref4 : r > ref4; j = 0 <= ref4 ? ++r : --r) {
      scantwo_data.int[i][j] = scantwo_data.full[i][j] - scantwo_data.add[i][j];
    }
  }
  // delete the non-symmetric versions
  scantwo_data.lod = null;
  scantwo_data.lodv1 = null;
  scantwo_data.max = {};
  ref5 = ["full", "add", "fv1", "av1", "int"];
  for (s = 0, len = ref5.length; s < len; s++) {
    i = ref5[s];
    scantwo_data.max[i] = d3panels.matrixMax(scantwo_data[i]);
  }
  return scantwo_data;
};

lod_for_heatmap = function lod_for_heatmap(scantwo_data, left, right) {
  var i, j, l, o, ref, ref1, thelod, z;
  // make copy of lod
  z = scantwo_data.full.map(function (d) {
    return d.map(function (dd) {
      return dd;
    });
  });
  for (i = l = 0, ref = z.length; 0 <= ref ? l < ref : l > ref; i = 0 <= ref ? ++l : --l) {
    for (j = o = 0, ref1 = z.length; 0 <= ref1 ? o < ref1 : o > ref1; j = 0 <= ref1 ? ++o : --o) {
      thelod = j < i ? right : left;
      z[i][j] = scantwo_data[thelod][i][j] / scantwo_data.max[thelod] * scantwo_data.max["full"];
    }
  }
  return z; // return the matrix we created
};