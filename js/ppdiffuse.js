/* Given a polypeptide input, calculates its voltage-dependent translocation
probability and dynamics for output in a webpage or downloadable format.
This is done by integrating the Smoluchowski equation over an arbitrary potential
with voltage-dependent electrostatic forces, length-dependent entropic forces, and
a barrier to translocation that can be used for molecules that bind to membranes
before presentation to the nanopore. */

var asynseqN2C = 'MDVFMKGLSK AKEGVVAAAE KTKQGVAEAA GKTKEGVLYV GSKTKEGVVH GVATVAEKTK EQVTNVGGAV VTGVTAVAQK TVEGAGSIAA ATGFVKKDQL GKNEEGAPQE GILEDMPVDP DNEAYEMPSE EGYQDYEPEA'
//var asynseqN2C = 'QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQEQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ'
var asynseq = asynseqN2C.split('').reverse().join('');
var kT = 1.38e-23*300*1e12*1e9; // units are pN nm
var e = 1.6e-19*1e12*1e9*1e-3; // sets voltage units in mV
//var laa = 0.4; // length per amino acid in nm; now an input
//var b = 0.6; // Kuhn length in nm; now an input

function cumsum(array, dx) {
	//cumulative sum of an array. Like an integral.
	
	var dx = dx || 1.0
	
	array = array.reduce(function(r, a) {
		if (r.length > 0) a += r[r.length - 1];
		r.push(a);
        return r;
		}, []);
	
	/*runsum = 0;
	for (var i = 0; i<array.length; i++) {
		runsum+=array[i];
		array[i] = runsum;
	}*/
	
	return math.multiply(array, dx)
}

function gaussian(x, x0, sigma) {
	// return a normalized gaussian	with the appropriate x0 and sigma
	return math.multiply(math.pow(2*math.PI*math.pow(sigma,2.0), -0.5), gaussian_unnorm(x, x0, sigma));
}

function gaussian_unnorm(x, x0, sigma) {
	// return an unnormalized gaussian (value = 1 at x-x0 = 0) with the appropriate x0 and sigma
	return math.exp(math.divide(math.subtract(0, math.dotPow(math.subtract(x, x0),2.0)), 2.0*math.pow(sigma, 2.0)));
}

function erf(x, x0, sigma, dx) {
	var dx = dx || 1.0
	return cumsum(gaussian(x, x0, sigma), dx)
}

function US(x, laa) {
	
	maxx = x[x.length-1]
	maxxp = maxx + 2*laa
	pre = 0.59*kT
	xplaa = math.add(x, laa)
	x2 = math.divide(xplaa, maxxp)
	
	//return math.multiply(pre, math.divide(math.subtract(1, math.multiply(x2, 2)), math.multiply(x2, math.subtract(1, x2)))) this is the forces
	return math.multiply(pre, math.add(math.log(x2), math.log(math.subtract(1, x2))))
	
}

function US1(x, Lp, Ltether, b) {

	var dx = (x[x.length-1] - x[0])/(x.length-1);
	var delta = 1e-3;
	
	function isPositive (v) {
		return v > 0;
	}
	
	function replacemin (a) {
		var mina = math.min(math.filter(a, isPositive))
		return math.map(a, function(v) {return (isPositive(v)) ? v : mina})
	}

	var L = x[x.length-1];
	var nt = math.divide(math.subtract(L-Lp/2, x), b);
	var nc = math.divide(math.subtract(x, Lp/2), b);
	var crit = math.and(math.map(nt, isPositive), math.map(nc, isPositive));
    nc = replacemin(nc);
	nt = replacemin(nt);
	
	//var Unt = math.erf(math.multiply(delta/b, math.sqrt(math.divide(3/2, nt))))
	var Unt = math.erf(math.multiply(delta/b, math.sqrt(math.multiply(3/2, math.dotPow(nt, -1)))))
	var ncpre = math.dotPow(math.multiply(2/3, math.pi, b*b, nc), -3/2)
	var ncpreexp = math.dotPow(math.multiply(-2/3, b*b, nc), -1)
	var Unc = math.dotMultiply(ncpre, math.subtract(math.exp(math.multiply(ncpreexp, math.pow(Ltether-delta, 2))), math.exp(math.dotMultiply(ncpreexp, math.pow(Ltether+delta, 2)))))

	var Un = math.dotMultiply(Unc, Unt)
	var minUng0 = math.min(math.filter(Un, isPositive))
	Un = math.map(Un,  function(v) {return (v<=0) ? minUng0 : v})
	
	return [math.multiply(-1., math.log(Un.reverse())), crit]
	
	
}

function seq2charge(seq, N2C) {
	// converts string seq into array of charge values
	// N2C is a boolean that determines which end is the N terminus;
	// true if residues are N->C, false if C->N.
	var numarray = seq.split('')		// convert to string
	N2C = N2C || false
	//match to amino acid type
	for (var i=0; i<numarray.length; i++) {
		switch (numarray[i]) {
			case 'X':
				numarray[i] = -2.0;
				break;
			case 'D':
			case 'E':
				numarray[i] = -1.0;
				break;
			case 'K':
			case 'R':
			case 'H':
				numarray[i] = 1.0;
				break;
			default:
				numarray[i] = 0.0;
		}
	}
	switch (N2C) {
		case true:
			numarray[0] = numarray[0] + 1;
			numarray[numarray.length-1] = numarray[numarray.length-1] - 1;
			break;
		default:
			numarray[0] = numarray[0] - 1;
			numarray[numarray.length-1] = numarray[numarray.length-1] + 1;
	}
			
	return numarray;
}

function seq2chargex(seq, N2C, x, laa, sigma) {
	// converts string seq into array of charge values
	// N2C is a boolean that determines which end is the N terminus;
	// true if residues are N->C, false if C->N.
	
	// x is a uniformly spaced range of points onto which seq is mapped
	// sigma is the width of the Gaussian convolution function
	var numarray = seq2charge(seq, N2C);
	var dx = (x[x.length-1] - x[0])/(x.length-1);
	//var dx = (x[x.length-1] - x[0])/numarray.length;
	y = math.zeros(x.length).toArray();
	for (var i = 0; i<numarray.length; i++) {
		dy = gaussian(x, (i + 0.5)*laa, sigma);
		y = math.add(y, math.multiply(numarray[i], dy, 1/(math.sum(dy)*dx)))
		//document.write(x[i], ' ', numarray[i], ' ', y[i], '<br>')
	}
	//document.write(y)
	return y
}

function applyEOF(chargedensity, eofm, eofb) {
	return math.map(chargedensity, function (c) {return eofm*c + eofb})
}


function interp(xf, xi, yi) {
	// 1D linear interpolation for a single point
	yf = [];
	if (xf.length == undefined) {xf = [xf];}
	for (var i=0; i<xf.length; i++) {
		xl = math.filter(xi, (function lt(x) {return x<=xf[i]})).slice(-1)[0];
		yl = yi[xi.indexOf(xl)];
		if (xl == xf[i]) {
			yf.push(yl);
		} else {
			xr = math.filter(xi, (function gt(x) {return x>xf[i]}))[0];
			yr = yi[xi.indexOf(xr)];
			yf.push(yl + (yr-yl)*(xf[i]-xl)/(xr-xl));
		}
	}
	return (xf.length==1) ? yf[0] : yf
}

function barrier(V, x, chargedensity, laa) {

	//var protseq = seq2charge(seq, N2C);
	//var x = math.multiply(math.range(0, protseq.length).toArray(), laa);
	var dx = (x[x.length-1] - x[0])/(x.length-1);	
	Ue = math.multiply(cumsum(chargedensity, dx), e*V/kT, 1.);
	//document.write(math.min(Ue), '<br>')
	//U = math.add(Ue, math.multiply(1., US(x, laa)), math.multiply(1*Eb, erf(x, xb, wb, dx)), math.multiply(1*fu, x), math.multiply(25, gaussian(x, 5., wb/2.3)))
	U = math.add(Ue, math.multiply(1., US(x, laa)), get_barrier_components(x, dx))
	//U = math.add(Ue, math.multiply(1., US1(x, 0.6, 4.0, 2.0)), math.multiply(0*Eb, erf(x, xb, wb, dx)), math.multiply(1*fu, x))
	U = math.subtract(U, U[0])
	return U
}

function barrier1(V, x, chargedensity, laa, b, Lp, Ltether) {

	//var protseq = seq2charge(seq, N2C);
	//var x = math.multiply(math.range(0, protseq.length).toArray(), laa);
	var dx = (x[x.length-1] - x[0])/(x.length-1);	
	Ue = math.multiply(cumsum(chargedensity, dx), e*V/kT, 1.);
	//document.write(math.min(Ue), '<br>')	
	//U = math.add(Ue, math.multiply(1., US(x, laa)), math.multiply(1*Eb, erf(x, xb, wb, dx)), math.multiply(1*fu, x))
	entropy = US1(x, Lp, Ltether, b)
	//U = math.add(Ue, math.multiply(1., entropy[0]), math.multiply(1*Eb, erf(x, xb, wb, dx)), math.multiply(1*fu, x))
	U = math.add(Ue, math.multiply(1., entropy[0]), get_barrier_components(x, dx))
	x = math.filter(x, function (v, i) {return entropy[1][i]})
	x = math.subtract(x, x[0])
	U = math.filter(U, function (v, i) {return entropy[1][i]})
	U = math.subtract(U, U[0])
	return [x, U]
}

function get_barrier_components(x, dx) {
	U = math.multiply(0., x)
	$('[id^=panel-]').each(function() {
		panelid = $(this).attr('id')
		paneltype = $(this).attr('paneltype')
		U = math.add(U, elements[$(this).attr('paneltype')].U(panelid, x, dx))
		})
	return U
}

function calcEscapeTime1x(x, x0, U, D, dx) {

	//everything is reversed
	//console.log(x0)
	x0 = x[x.length-1] - x0;
	var xg = math.filter(x, (function gt(x) {return x>=x0}));
	var idx = x.indexOf(xg[0]);
	var ridx = math.index(math.range(idx, x.length, 1, false))
	//console.log(x0, x, ridx.toArray())
		
	Ur = U.slice().reverse() //required because calculation assumes reflecting boundary on left, while it's actually on the right
	
	var innerintegral = cumsum(math.exp(math.multiply(-1., Ur)), dx)
	var outerintegral = math.multiply(math.sum(math.dotMultiply(math.exp(math.subset(Ur, ridx)), math.subset(innerintegral, ridx))), dx)
	
	return math.divide(outerintegral, D)
}

function calcEscapeTime2(x, x0, U, D, laa) {
	// calculate passage times
	//OBSOLETE
//	U = math.subtract(U, U[0]);
	Ur = U.slice().reverse(); // slice is necessary not to reverse in place

	var PiR = math.divide(cumsum(math.exp(U), laa), (math.sum(math.exp(U)) * laa));
	var thetaRiterm = math.multiply(math.dotMultiply(math.exp(U), cumsum(math.dotMultiply(PiR, math.exp(math.multiply(U, -1))), laa)), 1/D);
	var thetaRterm = math.multiply((math.sum(thetaRiterm)*laa),cumsum(math.exp(U), laa), 1/(math.sum(math.exp(U))*laa));
	var thetaR = math.subtract(thetaRterm,cumsum(thetaRiterm, laa));

	var PiL = math.divide(cumsum(math.exp(Ur), laa), (math.sum(math.exp(Ur))*laa)).reverse();
	var thetaLiterm = math.multiply(1/D, math.dotMultiply(math.exp(Ur), cumsum(math.dotMultiply(PiL.slice().reverse(), math.exp(math.multiply(Ur, -1))), dx)));
	var thetaLterm = math.multiply((math.sum(thetaLiterm)*laa), cumsum(math.exp(Ur), laa), 1/(math.sum(math.exp(Ur))*laa));
	var thetaL = math.subtract(thetaLterm, cumsum(thetaLiterm, laa)).reverse();
	
	var pisideL = interp(x0, x, PiL);
	var pisideR = interp(x0, x, PiR);
	
	var tescL = interp(x0, x, math.dotMultiply(thetaL, math.dotPow(PiL, -1)));
	var tescR = interp(x0, x, math.dotMultiply(thetaR, math.dotPow(PiR, -1)));
	
	var tesc = interp(x0, x, math.add(thetaL, thetaR));
	
//	for (var i=0; i<x.length; i++) {document.write(i, '\t', x[i], '\t', math.dotMultiply(thetaL, math.dotPow(PiL, -1))[i], '<br>');}
	
	return [pisideL, pisideR, tescL, tescR, tesc]


}

function calcEscapeTime2x(x, x0, U, D, dx) {
	// calculate passage times
//	U = math.subtract(U, U[0]);
	Ur = U.slice().reverse(); // slice is necessary not to reverse in place

	var PiR = math.divide(cumsum(math.exp(U), dx), (math.sum(math.exp(U)) * dx));
	var thetaRiterm = math.multiply(math.dotMultiply(math.exp(U), cumsum(math.dotMultiply(PiR, math.exp(math.multiply(U, -1))), dx)), 1/D);
	var thetaRterm = math.multiply((math.sum(thetaRiterm)*dx),cumsum(math.exp(U), dx), 1/(math.sum(math.exp(U))*dx));
	var thetaR = math.subtract(thetaRterm,cumsum(thetaRiterm, dx));

	var PiL = math.divide(cumsum(math.exp(Ur), dx), (math.sum(math.exp(Ur))*dx)).reverse();
	var thetaLiterm = math.multiply(1/D, math.dotMultiply(math.exp(Ur), cumsum(math.dotMultiply(PiL.slice().reverse(), math.exp(math.multiply(Ur, -1))), dx)));
	var thetaLterm = math.multiply((math.sum(thetaLiterm)*dx), cumsum(math.exp(Ur), dx), 1/(math.sum(math.exp(Ur))*dx));
	var thetaL = math.subtract(thetaLterm, cumsum(thetaLiterm, dx)).reverse();
	
	var pisideL = interp(x0, x, PiL);
	var pisideR = interp(x0, x, PiR);
	
	var tescL = interp(x0, x, math.dotMultiply(thetaL, math.dotPow(PiL, -1)));
	var tescR = interp(x0, x, math.dotMultiply(thetaR, math.dotPow(PiR, -1)));
	
	var tesc = interp(x0, x, math.add(thetaL, thetaR));
	
//	for (var i=0; i<x.length; i++) {document.write(i, '\t', x[i], '\t', math.dotMultiply(thetaL, math.dotPow(PiL, -1))[i], '<br>');}
	
	return [pisideL, pisideR, tescL, tescR, tesc]


}

function minInRange(x, y, xl, xr) {
	// finds position of function (x, y) minimum between xl and xr
	// requires at least 3 total points
	
	idxs = math.range(0, x.length).toArray()
	xcrit = math.filter(x, (function crit(x) {return ((x>=xl) && (x<=xr))}))
    idxrange = math.range(x.indexOf(xcrit[0]), x.indexOf(xcrit[xcrit.length-1])+1)
	xsub = math.subset(x, math.index(idxrange))
	ysub = math.subset(y, math.index(idxrange))
	
	//if result is a single point, return it; otherwise find minimum value
	return (xsub.length==undefined) ? [xsub, ysub] : [xsub[ysub.indexOf(math.min(ysub))], math.min(ysub)]
}

function makeChargeDensityControls(target_id) {

	var exportControls = d3.select("#" + target_id).append('div')
			.classed("exportControls", true)
			
		exportControls
		.append("button")
			.text("Export")
			.classed("ui-button", true)
			.attr("id", "btnExportChargeDensity")
			.on("click", function () {
					//format data to saveData
					var output = formatOutput(chart3);
					saveData(d3.tsvFormatRows(output), "chargedensity.txt", "text/tab-separated-values");
				})	
}

function makePotentialControls(target_id_right, target_id_left) {

	var exportControls = d3.select("#" + target_id_right).append('div')
			.classed("exportControls", true)
			
		exportControls
		.append("button")
			.text("Export")
			.classed("ui-button", true)
			.attr("id", "btnExportPotential")
			.on("click", function () {
					//format data to saveData
					var output = formatOutput(chart2);
					//console.log(output, d3.tsvFormatRows(output))
					saveData(d3.tsvFormatRows(output), "potential.txt", "text/tab-separated-values");
				})
	var injectControls = d3.select("#" + target_id_left).append('div')
				.classed("injectControls", true)
				
		injectControls
			.append("input")
			.attr("id", "showInjectionPoints")
			.attr("type", "checkbox")
			//.style("zIndex", 0)
			.property("checked", false)
			.on("change", function() {
				update_Uplot();
				fitPlots();
			})		
		injectControls
			.append("label")
			.text("Show injection points")
			
}

function formatOutput(targetchart) {
	var labels = targetchart.options().series;
	var data = targetchart.source_data();
	var row = [];
	var headerrow = ["x (nm)"];
	var output = [];
	for (var i = 0; i < labels.length; i++) {
		headerrow.push(labels[i]["label"])
	}
	output.push(headerrow)
	for (var j = 0; j < data[0].length; j++) {
		row.push(data[0][j][0]);
		for (var i = 0; i < data.length; i++) {
			row.push(data[i][j][1]);
		}
		output.push(row);
		row = [];
	}
	return output;
}

function makeGeneralProfileControls(target_id) {
	
	  changefunc = update_plots;
	
	  var genControls = d3.select("#" + target_id).append('div')
	    .classed("gencontrols", true)
      genControls
//		.append("label")
//          .text("General parameters: ")
//		  .style("font-weight", "bold")
		.append("li")
		  .text("Diffusion constant (\u03BCm\u00B2/s): ")
		  .style("font-weight", "normal")
        .append("input")
		  .attr("id", "D")
		  .attr("save", true)
          .attr("type", "text")
          .attr("name", "D")
          .attr("value", 0.4)
		  .on("change", changefunc)
/*	  genControls
		.append("li")
		  .text(" constant offset force (pN): ")
		.append("input")
		  .attr("id", "fu")
		  .attr("save", true)
          .attr("type", "text")
          .attr("name", "fu")
          .attr("value", 0.0)
		  .on("change", changefunc) */
}

function makePoreControls(target_id) {
	
	 changefunc = update_plots;
	
	 var poreControls = d3.select("#" + target_id).append('div')
	    .classed("porecontrols", true)
	  poreControls
	  	/*.append("label")
          .text("Pore parameters: ")
		  .style("font-weight", "bold")*/
		.append("li")
		  .text("Length (nm): ")
		  .style("font-weight", "normal")
		.append("input")
		  .attr("id", "Lp")
		  .attr("save", true)
          .attr("type", "text")
          .attr("name", "Lp")
          .attr("value", 1.2)
		  .on("change", changefunc)
	  poreControls
	  	  .append("li")
		  .text("Electroosmotic slope: ")
	      .append("input")
		  .attr("id", "EOFm")
		  .attr("save", true)
          .attr("type", "text")
          .attr("name", "EOFm")
          .attr("value", 0.908)
		  .on("change", changefunc)
	  poreControls
   	    .append("li")
		  .text("Electroosmotic intercept (e/nm): ")
		.append("input")
		  .attr("id", "EOFb")
		  .attr("save", true)
          .attr("type", "text")
          .attr("name", "EOFb")
          .attr("value", -0.18)
		  .on("change", changefunc)
      
}

function makeSequenceControls(target_id, initseq) {
	
	changefunc = update_plots;
	
	var seqControls = d3.select("#" + target_id).append('div')
		.classed("seqcontrols", true)
	  
	  seqControls
        .append("label")
		  .attr("id", "seqlabel")
		  .text("Sequence: ")
		.append("button")
			.text("reverse")
			.classed("ui-button ui-widget", true)		
			.attr("id", "reversebutton")
			.on("click", function () {
				$("#sequence").val($("#sequence").val().split('').reverse().join(''));
				($("#seqdirN2C").prop("checked")) ? $("#seqdirC2N").prop("checked", true) : $("#seqdirN2C").prop("checked", true)
				changefunc();
				})
	  seqControls.append("br")
	  seqControls
        .append("textarea")
		  .attr("id", "sequence")
		  .attr("save", true)
          .attr("type", "text")
          .attr("name", "sequence")
          .attr("text", initseq)
		  //.val(initseq)
		  .on("change", changefunc)
		  $("#sequence").val(initseq)

		seqControls.append("br")
		seqControls
        .append("label")
		.attr("id", "seqlabel2")
        .text(" Direction:")
        .append("input")
          .attr("type", "radio")
		  .attr("save", true)
		  .attr("name", "seqdir")
		  .attr("id", "seqdirN2C")
          .property("checked", false)
          .attr("value", "N2C")
          .on("change", changefunc)
		seqControls
        .append("label")
		.attr("id", "seqlabel3")
        .text("N->C")
        .append("input")
          .attr("type", "radio")
		  .attr("save", true)
          .property("checked", true)
		  .attr("name", "seqdir")
		  .attr("id", "seqdirC2N")		  
          .attr("value", "C2N")
          .on("change", changefunc)
		 seqControls
		 .append("label")
		 .attr("id", "seqlabel4")
		 .text("C->N") 
		 seqControls2 = seqControls.append("div")
		 seqControls2
		  .append("label")
		  .text("Length per amino acid (nm): ")
	      .append("input")
			.attr("id", "laa")
			.attr("save", true)
			.attr("type", "text")
			.attr("name", "laa")
			.attr("value", 0.4)
			.on("change", changefunc)
		  seqControls2
		  .append("label")
		  .text(" Kuhn length (nm): ")
	      .append("input")
			.attr("id", "b_kuhn")
			.attr("save", true)
			.attr("type", "text")
			.attr("name", "b_kuhn")
			.attr("value", 0.6)
			.on("change", changefunc)
}

function makeGeneralAccordionControls(target_id) {

	changefunc = update_plots;
	
	var accControls = d3.select("#" + target_id + "_topleft")
	
	accControls
		.append("button")
			.text("Add new: ")
			.classed("ui-button", true)		
			.attr("id", "addelement")
			.style("background", "#ffff99")
			.on("click", function () {
				
				var newlist = new Array()
				$('.ui-accordion-header').each(function() {newlist.push($(this).attr('id').split('-').slice(-1)[0])})
				if (newlist.length == 0) {newlist = ["-1"]}
				var nextnum = (Number(newlist.slice(-1)[0]) + 1).toString()
				var panelid = "panel-" + nextnum
				var newtype = $("#new_element_type").val()
	
				makeNewPanel(panelid, newtype)
				
				})
	accControls
		/*.append("label")
          .text(" Type: ")
		  //.style("font-weight", "bold")*/
		.append("select")
		  .classed("ui-button", true)
		  //.style("background", "#ffff99")
		  .attr("id", "new_element_type")
		  .selectAll("option")
		  .data($.map(elements, function (v, k) {return v.label;}))
		  .enter().append("option")
		    .text(function (d) {return d;})
			
		//$("#new_element_type").selectmenu()

	var accControls = d3.select("#" + target_id + "_topright")
	
	accControls
		.append("label")
		.text("Configuration: ")
		.append("button")
		.text("Save")
		.classed("ui-button", true)		
		.on("click", saveConfig)
	accControls
		.append("input")
		.attr("style", "display:none;")
		.attr("type", "file")
		.attr("id", "loadconfig")
		.text("Load")
		.on("change", function() {
			loadConfig();
			$(this).val("");
		})
	accControls
		.append("button")
		.text("Load")
		.classed("ui-button", true)
		.on("click", function() {$("#loadconfig").click()})
	accControls.append("p")

}

function genPanelControls(panelid, fields, changefunc) {
			var ctrl = d3.select("#" + panelid)
			$.each(fields, function(key, field) {
				ctrl.append("li")
				  .text(field.label + ": ")
				  //.style("font-weight", "normal")
				.append("input")
				  .attr("id", key + "-" + panelid)
				  .attr("save", true)
				  .attr("type", "text")
				  .attr("name", key)
				  .attr("value", field.defaultValue)
				  .on("change", changefunc)
			})
			
}

function makeNewPanel(panelid, paneltype) {
	
	changefunc = update_plots;
	
	nextnum = panelid.split('-').slice(-1)[0]
	
	for (var key in elements) {
		if (elements[key].label==paneltype) {
			newtype = key;
			break;
		}
	}
	
	var panelContent = "<h3>" + paneltype + "</h3><div paneltype=" + newtype + " id=" + panelid + "></div>";
	$('#accordion').append(panelContent)

	elements[newtype].createPanelControls(panelid, changefunc);
	
	newControls = d3.select("#" + panelid)
		.append('div')
		.append("button")
			.text("Delete")
			.classed("ui-button ui-widget", true)		
			.attr("id", "removeelement")
			.on("click", function () {

				var parent = $(this).parentsUntil("#accordion").last();
				parent.prev().remove();
				parent.remove();
				changefunc();
				
				})
	$("#accordion").accordion('destroy').accordion({collapsible: true, active: -1});
	changefunc();
}

function makeEntropyControls(target_id) {

	  changefunc = update_plots;

      var profControls = d3.select("#" + target_id).append('div')
          .classed("profcontrols", true)
      
	  /*profControls.append("button")
        .text("Calculate")
        .classed("ui-button ui-corner-all ui-widget", true)		
        .on("click", update_plots)*/
      
      var modeControls = d3.select("#" + target_id).append('div')
         .classed("modecontrols", true)
		 
      modeControls
		.append("input")
		.attr("id", "tethercheck")
		.attr("save", true)
		.attr("type", "checkbox")
		.property("checked", false)
		.on("change", function () {
			$("input[name*='tethergroup']").prop("disabled", !($("#tethercheck").prop("checked")));
			var oldval = $("#Yplotquantity").val()
			labels = ($("#tethercheck").prop("checked")) ? tlabelmap : labelmap
			d3.select("#Yplotquantity").selectAll("option")
  			  .remove()
			d3.select("#Yplotquantity").selectAll("option")
			  .data($.map(labels, function (v, k) {return v;}))
			  .enter().append("option")
			    .text(function (d) {return d;})

			var newval = labels["tesc"]
			for (var key in labels) {
				if (labels[key]==oldval) {
					newval = oldval;
					break;
				}
			}				
			$("#Yplotquantity").val(newval)
			changefunc();
		})
	  modeControls
	   .append("label")
		.text("Polypeptide tethered? ")
	  modeControls
		.append("label")
		  .text("Tether length (nm): ")
		  .attr("name", "tethergroup0")
		.append("input")
		  .attr("id", "Ltether")
		  .attr("save", true)
		  .attr("type", "text")
		  .attr("name", "tethergroup1")
		  .attr("value", 2.0)
		  .property("disabled", true)
		  .on("change", changefunc)

}

function makeVoltageControls(target_id) {

	changefunc = update_plots;

	var voltageControls = d3.select("#" + target_id).append('div')
	    .classed("gencontrols", true)
     voltageControls
		.append("label")
          .text("Voltage range (mV): ")
		  //.style("font-weight", "bold")
		.append("label")
		  .text("min: ")
		  .style("font-weight", "normal")
        .append("input")
		  .attr("id", "minV")
		  .attr("save", true)
          .attr("type", "text")
          .attr("name", "minV")
          .attr("value", 10.)
		  .on("change", changefunc)
	  voltageControls
		.append("label")
		  .text(" max: ")
		.append("input")
		  .attr("id", "maxV")
		  .attr("save", true)
          .attr("type", "text")
          .attr("name", "maxV")
          .attr("value", 60.)
		  .on("change", changefunc)
	  voltageControls
		.append("label")
		  .text(" step: ")
		.append("input")
		  .attr("id", "stepV")
		  .attr("save", true)
          .attr("type", "text")
          .attr("name", "stepV")
          .attr("value", 2.5)
		  .on("change", changefunc)
		 
	  var minControls = d3.select("#" + target_id).append('div')
	    .classed("mincontrols", true)
      minControls
		.append("label")
		  .text("Injection point (nm):")
		  //.style("color", optsx.color1)
	  minControls
		.append("label")
		  .text(" at potential minimum between ")
        .append("input")
		  .attr("id", "x0min")
		  .attr("save", true)
          .attr("type", "text")
          .attr("name", "x0min")
          .attr("value", 5.)
		  .on("change", changefunc)
	  minControls.append("label")
          .text(" and ")		  
		.append("input")
		  .attr("id", "x0max")
		  .attr("save", true)
          .attr("type", "text")
          .attr("name", "x0max")
          .attr("value", 20.)
		  .on("change", changefunc)
		  
}

function makeVplotControls(target_id_topleft, target_id_topcenter, target_id_topright, target_id_left, target_id_right, labels) {
	
	var changefunc = function () {
		update_Vplot();
		fitPlots();
	}
	
	var scaleControls = d3.select("#" + target_id_left).append('div')
	    .classed("scalecontrols", true)
	  /*scaleControls
		.append("label")
		  .text(" Y axis: ")
		.append("select")
		  .data([)
		  .enter().append("label")
		  .text(function (d) {return labels[d]})
		  .attr("id", "Yid")
          .attr("type", "option")
          .attr("name", "Yid")
          .attr("value", "tesc")
		  .on("change", changefunc)
		  */
	  scaleControls
		.append("label")
          .text("Y axis: ")
		  //.style("font-weight", "bold")
		.append("select")
		  .attr("id", "Yplotquantity")
		  .on("change", changefunc)
		  .selectAll("option")
		  .data($.map(labels, function (v, k) {return v;}))
		  .enter().append("option")
		    .text(function (d) {return d;})
	  	        
	   scaleControls
		.append("select")
		  .attr("id", "Yplotscale")
		  .on("change", changefunc)
		  .selectAll("option")
		  .data(["linear", "log"])
		  .enter().append("option")
		    .text(function (d) {return d;})
		
		$("#Yplotquantity").val(labels["tesc"]);
		$("#Yplotscale").val("log");
		
	var fileControls = d3.select("#" + target_id_topleft)
		
		fileControls
		.append("label")
			.text("Load data file: ")
		.append("input")
//			.classed("ui-button", true)
			.attr("type", "file")
			.on("change", loadData)
			.attr("id","datafile")
			
	var fileControls2 = d3.select("#" + target_id_topright)
			
		fileControls2
		.append("button")
			//.classed("ui-button", true)
			.text("Clear data")
			.attr("id", "btnClearData")
			.on("click", function () {
				$("#datafile").val("");
				loadData();
				})
	
	var fileControls3 = d3.select("#" + target_id_topcenter)
			
		fileControls3
		.append("button")
			//.classed("ui-button", true)
			.text("Reverse polarity")
			.attr("id", "btnReverseV")
			.property("disabled", true)
			.on("click", function () {
				for (var i = 0; i < expdata.length; i++) {
					expdata[i][0] = -expdata[i][0];
					expdata[i][2]["xlower"] = -expdata[i][2]["xlower"];
					expdata[i][2]["xupper"] = -expdata[i][2]["xupper"];
				}
				changefunc();
				})
				
	var exportControls = d3.select("#" + target_id_right).append('div')
			.classed("exportControls", true)
			
		exportControls
		.append("button")
			.text("Export calculation")
			.classed("ui-button", true)
			.attr("id", "btnExportCalc")
			.on("click", function () {
					//format data to saveData
					var output = d3.tsvFormatRows(chart.source_data()[0]);
					saveData(output, "calculation.txt", "text/tab-separated-values");
				})
}

function loadData() {

	var changefunc = function () {
		update_Vplot();
		fitPlots();
	}

	var file = document.getElementById("datafile").files[0]; // only one file allowed
	if (file) {
		datafilename = file.name;
		var result = null;
		var reader = new FileReader();
		reader.onload = function(e) {

			r = d3.tsvParseRows(e.target.result, function(d) {
				return d.map(Number);
			})
			expdata = [];
			for (var i = 0; i<r.length; i++) {
				expdata.push([r[i][0], r[i][1], {"xlower": r[i][0], "xupper": r[i][0], "ylower": r[i][1]-r[i][2], "yupper": r[i][1] + r[i][2]}])
			}
			
			$("#btnReverseV").prop("disabled", false)
			//console.log(expdata);
			//Format for plotting
			changefunc();
		}
		reader.readAsText(file);
	} else {
		expdata = [];
		$("#btnReverseV").prop("disabled", true)
		changefunc();
	}
}

function saveData(data, fileName, type) {
	var type = type || "application/python";
	var blob = new Blob([data], {type: type}),
		url = window.URL.createObjectURL(blob);
	if (window.navigator.msSaveOrOpenBlob) { 
	  window.navigator.msSaveOrOpenBlob(blob, fileName); 
	} else {
	  var a = d3.select("body").append("a")
  //.style("display", "none")
	   .style("display", "none")
	   .attr("id", "savedata")
	   .attr("href", url)
	   .attr("download", fileName)
	   .node().click();
	  setTimeout(function() { window.URL.revokeObjectURL(url) }, 1000);
	}
}

function saveConfig() {
	// saves current configuration as a JSON file that can be reloaded later
	
	var panels = new Object();
	var data = new Object();

	// save each custom panel (id = panel-*)
	$('[id^=panel-]').each(function() {
		panels[$(this).attr('id')] = $(this).attr('paneltype');
	})

	// save each control marked with a "save" tag
	$("[save=true]").each(function() {
		  switch ($(this).attr('type')) {
			  case 'text':
				 data[$(this).attr('id')] = $(this).val();
				 break;
			  case 'radio':
			  case 'checkbox':
				 data[$(this).attr('id')] = $(this).prop('checked');
				 break;
			  default:
			}
			})
	
	output = JSON.stringify({panels: panels, data: data, expdata: expdata}, null, 1)
	
	saveData(output, "config.json", "application/json");
}

function loadConfig() {
	
	var file = $("#loadconfig")[0].files[0]; // only one file allowed
	if (file) {
		datafilename = file.name;
		var result = null;
		var reader = new FileReader();
		reader.onload = function(e) {
		
			// first parse JSON back into config, then call this function
			config = JSON.parse(e.target.result);

			// clear accordion of custom panels
			$("[id^=panel-]").each(function () {
				$(this).prev().remove();
				$(this).remove();
			})
			$("#accordion").accordion('destroy').accordion({collapsible: true, active: false});
			
			// recreate custom panels
			for (var key in config["panels"]) {
				var paneltype = elements[config["panels"][key]].label
				makeNewPanel(key, paneltype)
			}
			
			// fill in all the fields
			for (var key in config["data"]) {
				switch ($("#" + key).attr('type')) {
					case "text":
						$("#" + key).val(config["data"][key])
						break;
					case "radio":
					case "checkbox":
						$("#" + key).prop("checked", config["data"][key])
						break;
					default:
				}
			}

			// read in experiment data if saved ("undefined" check is for backwards compatibility)
			expdata = (config["expdata"] == undefined) ? [] : config["expdata"];

			update_plots()
		}
		reader.readAsText(file);
	}
}

