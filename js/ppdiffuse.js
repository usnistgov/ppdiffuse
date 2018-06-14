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
var laa = 0.4; // length per amino acid in nm
var b = 0.6; // Kuhn length in nm

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
	return math.multiply(math.pow(2*math.PI*math.pow(sigma,2.0), -0.5),math.exp(math.divide(math.subtract(0, math.dotPow(math.subtract(x, x0),2.0)), 2.0*math.pow(sigma, 2.0))));
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

function US1(x, Lp, Ltether) {

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

function barrier(V, x, chargedensity, laa, Eb, xb, wb, fu) {

	//var protseq = seq2charge(seq, N2C);
	//var x = math.multiply(math.range(0, protseq.length).toArray(), laa);
	var dx = (x[x.length-1] - x[0])/(x.length-1);	
	Ue = math.multiply(cumsum(chargedensity, dx), e*V/kT, 1.);
	//document.write(math.min(Ue), '<br>')
	U = math.add(Ue, math.multiply(1., US(x, laa)), math.multiply(1*Eb, erf(x, xb, wb, dx)), math.multiply(1*fu, x))
	//U = math.add(Ue, math.multiply(1., US1(x, 0.6, 4.0, 2.0)), math.multiply(0*Eb, erf(x, xb, wb, dx)), math.multiply(1*fu, x))
	U = math.subtract(U, U[0])
	return U
}

function barrier1(V, x, chargedensity, laa, Lp, Ltether, Eb, xb, wb, fu) {

	//var protseq = seq2charge(seq, N2C);
	//var x = math.multiply(math.range(0, protseq.length).toArray(), laa);
	var dx = (x[x.length-1] - x[0])/(x.length-1);	
	Ue = math.multiply(cumsum(chargedensity, dx), e*V/kT, 1.);
	//document.write(math.min(Ue), '<br>')	
	//U = math.add(Ue, math.multiply(1., US(x, laa)), math.multiply(1*Eb, erf(x, xb, wb, dx)), math.multiply(1*fu, x))
	entropy = US1(x, Lp, Ltether)
	U = math.add(Ue, math.multiply(1., entropy[0]), math.multiply(1*Eb, erf(x, xb, wb, dx)), math.multiply(1*fu, x))
	x = math.filter(x, function (v, i) {return entropy[1][i]})
	x = math.subtract(x, x[0])
	U = math.filter(U, function (v, i) {return entropy[1][i]})
	U = math.subtract(U, U[0])
	return [x, U]
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
	return (xsub.length==undefined) ? xsub : xsub[ysub.indexOf(math.min(ysub))]
}

function makeChargeDensityControls(target_id) {

	var exportControls = d3.select("#" + target_id).append('div')
			.classed("exportControls", true)
			
		exportControls
		.append("button")
			.text("Export")
			.attr("id", "btnExportChargeDensity")
			.on("click", function () {
					//format data to saveData
					var output = formatOutput(chart3);
					saveData(d3.tsvFormatRows(output), "chargedensity.txt", "text/tab-separated-values");
				})	
}

function makePotentialControls(target_id) {

	var exportControls = d3.select("#" + target_id).append('div')
			.classed("exportControls", true)
			
		exportControls
		.append("button")
			.text("Export")
			.attr("id", "btnExportPotential")
			.on("click", function () {
					//format data to saveData
					var output = formatOutput(chart2);
					//console.log(output, d3.tsvFormatRows(output))
					saveData(d3.tsvFormatRows(output), "potential.txt", "text/tab-separated-values");
				})	
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


function makeProfileControls(target_id, initseq) {

	  changefunc = update_plots;

      var profControls = d3.select("#" + target_id).append('div')
          .classed("profcontrols", true)
      
	  /*profControls.append("button")
        .text("Calculate")
        .classed("ui-button ui-corner-all ui-widget", true)		
        .on("click", update_plots)*/
      
	  var seqControls = profControls.append('div')
          .classed("seqcontrols", true)
	  
	  seqControls
        .append("label")
		  .attr("id", "seqlabel")
          .text("Sequence: ")
	  seqControls
        .append("input")
		  .attr("id", "sequence")
          .attr("type", "text")
          .attr("name", "sequence")
          .attr("value", initseq)
		  .on("change", changefunc)
	   seqControls.append("button")
        .text("reverse")
        .classed("ui-button ui-widget", true)		
		.attr("id", "reversebutton")
        .on("click", function () {
			$("#sequence").val($("#sequence").val().split('').reverse().join(''));
			($("#seqdirN2C").prop("checked")) ? $("#seqdirC2N").prop("checked", true) : $("#seqdirN2C").prop("checked", true)
			changefunc();
			})
			
	   seqControls
        .append("label")
		.attr("id", "seqlabel2")
        .text(" Direction:")
        .append("input")
          .attr("type", "radio")
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
          .property("checked", true)
		  .attr("name", "seqdir")
		  .attr("id", "seqdirC2N")		  
          .attr("value", "C2N")
          .on("change", changefunc)
		 seqControls
		 .append("label")
		 .attr("id", "seqlabel4")
		 .text("C->N")
		 
	  var genControls = d3.select("#" + target_id).append('div')
	    .classed("gencontrols", true)
      genControls
		.append("label")
          .text("General parameters: ")
		  .style("font-weight", "bold")
		.append("label")
		  .text("D (\u03BCm\u00B2/s): ")
		  .style("font-weight", "normal")
        .append("input")
		  .attr("id", "D")
          .attr("type", "text")
          .attr("name", "D")
          .attr("value", 0.4)
		  .on("change", changefunc)
	  genControls
		.append("label")
		  .text(" constant offset force (pN): ")
		.append("input")
		  .attr("id", "fu")
          .attr("type", "text")
          .attr("name", "fu")
          .attr("value", 0.0)
		  .on("change", changefunc)
	  
	  var poreControls = d3.select("#" + target_id).append('div')
	    .classed("porecontrols", true)
	  poreControls
	  	.append("label")
          .text("Pore parameters: ")
		  .style("font-weight", "bold")
		.append("label")
		  .text(" length (nm): ")
		  .style("font-weight", "normal")
		.append("input")
		  .attr("id", "Lp")
          .attr("type", "text")
          .attr("name", "Lp")
          .attr("value", 1.2)
		  .on("change", changefunc)
	  poreControls
	  	  .append("label")
		  .text(" Electroosmotic slope: ")
		.append("input")
		  .attr("id", "EOFm")
          .attr("type", "text")
          .attr("name", "EOFm")
          .attr("value", 0.908)
		  .on("change", changefunc)
		  
	  poreControls
	  	  .append("label")
		  .text(" intercept (e/nm): ")
		.append("input")
		  .attr("id", "EOFb")
          .attr("type", "text")
          .attr("name", "EOFb")
          .attr("value", -0.18)
		  .on("change", changefunc)
      
      var modeControls = d3.select("#" + target_id).append('div')
         .classed("modecontrols", true)
		 
      modeControls
		.append("input")
		.attr("id", "tethercheck")
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
		.append("label")
		  .text("Tether length (nm): ")
		  .attr("name", "tethergroup0")
		.append("input")
		  .attr("id", "Ltether")
		  .attr("type", "text")
		  .attr("name", "tethergroup1")
		  .attr("value", 2.0)
		  .property("disabled", true)
		  .on("change", changefunc)

	  
	  var barrierControls = d3.select("#" + target_id).append('div')
	    .classed("barriercontrols", true)
      barrierControls
		.append("label")
          .text("Barrier parameters: ")
		  .style("font-weight", "bold")
		.append("label")
		  .text("Barrier height (kT): ")
		  .style("font-weight", "normal")
        .append("input")
		  .attr("id", "Eb")
          .attr("type", "text")
          .attr("name", "Eb")
          .attr("value", 16.)
		  .on("change", changefunc)
	  barrierControls
		.append("label")
		  .text(" Barrier width (nm): ")
		.append("input")
		  .attr("id", "wb")
          .attr("type", "text")
          .attr("name", "wb")
          .attr("value", 3.5)
		  .on("change", changefunc)
	  barrierControls
		.append("label")
		  .text(" Barrier position (nm): ")
		.append("input")
		  .attr("id", "xb")
          .attr("type", "text")
          .attr("name", "xb")
          .attr("value", 12.)
		  .on("change", changefunc)
		  
	var voltageControls = d3.select("#" + target_id).append('div')
	    .classed("voltagecontrols", true)
      voltageControls
		.append("label")
          .text("Voltage range (mV): ")
		  .style("font-weight", "bold")
		.append("label")
		  .text("min: ")
		  .style("font-weight", "normal")
        .append("input")
		  .attr("id", "minV")
          .attr("type", "text")
          .attr("name", "minV")
          .attr("value", 10.)
		  .on("change", changefunc)
	  voltageControls
		.append("label")
		  .text(" max: ")
		.append("input")
		  .attr("id", "maxV")
          .attr("type", "text")
          .attr("name", "maxV")
          .attr("value", 60.)
		  .on("change", changefunc)
	  voltageControls
		.append("label")
		  .text(" step: ")
		.append("input")
		  .attr("id", "stepV")
          .attr("type", "text")
          .attr("name", "stepV")
          .attr("value", 2.5)
		  .on("change", changefunc)
		 
	  var minControls = d3.select("#" + target_id).append('div')
	    .classed("mincontrols", true)
      minControls
		.append("label")
          .text("Search for potential minimum between ")
        .append("input")
		  .attr("id", "x0min")
          .attr("type", "text")
          .attr("name", "x0min")
          .attr("value", 5.)
		  .on("change", changefunc)
	  minControls.append("label")
          .text(" nm and ")		  
		.append("input")
		  .attr("id", "x0max")
          .attr("type", "text")
          .attr("name", "x0max")
          .attr("value", 20.)
		  .on("change", changefunc)
	   minControls.append("label")
          .text(" nm")		  
	  
	  

		//$("#sequence").width($("#" + target_id).width());
//          .on("change", update_plots)
		
     /* 
      var fitControls = d3.select("#" + target_id).append('div')
        .classed("fit controls", true)
        .style("padding-top", "5px")
        
      fitControls.append("button")
        .text("start fit")

        .on("click", fit)
        
      fitControls.append("label")
        .text(" log:")
      
      fitControls.append("div")
        .append("pre")
        .classed("fit log", true)

      update_mode.call({value: "edit"});
      */
      //$(profControls.node()).controlgroup();
    }
	
function makeVplotControls(target_id_topleft, target_id_topright, target_id_left, target_id_right, labels) {
	
	var changefunc = update_Vplot
	
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
		
	var fileControls = d3.select("#" + target_id_topleft).append('div')
			.classed("fileControls", true)
		
		fileControls
		.append("label")
			.text("Load data file: ")
		.append("input")
//			.classed("ui-button ui-widget", true)
			.attr("type", "file")
			.on("change", loadData)
			.attr("id","datafile")
			
	var fileControls2 = d3.select("#" + target_id_topright).append('div')
			.classed("fileControls", true)
			
		fileControls2
		.append("button")
//			.classed("ui-button ui-widget", true)
			.text("Clear data")
			.attr("id", "btnClearData")
			.on("click", function () {
				$("#datafile").val("");
				loadData();
				})
				
	var exportControls = d3.select("#" + target_id_right).append('div')
			.classed("exportControls", true)
			
		exportControls
		.append("button")
			.text("Export calculation")
			.attr("id", "btnExportCalc")
			.on("click", function () {
					//format data to saveData
					var output = d3.tsvFormatRows(chart.source_data()[0]);
					saveData(output, "calculation.txt", "text/tab-separated-values");
				})
}

function loadData() {
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
			//console.log(expdata);
			//Format for plotting
			update_Vplot();
		}
		reader.readAsText(file);
	} else {
		expdata = [];
		update_Vplot();
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

