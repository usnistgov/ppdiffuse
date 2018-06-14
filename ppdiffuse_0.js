/* Given a polypeptide input, calculates its voltage-dependent translocation
probability and dynamics for output in a webpage or downloadable format.
This is done by integrating the Smoluchowski equation over an arbitrary potential
with voltage-dependent electrostatic forces, length-dependent entropic forces, and
a barrier to translocation that can be used for molecules that bind to membranes
before presentation to the nanopore. */

var asynseqN2C = 'MDVFMKGLSK AKEGVVAAAE KTKQGVAEAA GKTKEGVLYV GSKTKEGVVH GVATVAEKTK EQVTNVGGAV VTGVTAVAQK TVEGAGSIAA ATGFVKKDQL GKNEEGAPQE GILEDMPVDP DNEAYEMPSE EGYQDYEPEA'
var asynseq = asynseqN2C.split('').reverse().join('').replace(/\s/g, '');
var kT = 1.38e-23*300*1e12*1e9; // units are pN nm
var e = 1.6e-19*1e12*1e9*1e-3; // sets voltage units in mV
var laa = 0.4; // in nm

function cumsum(array, dx) {
	//cumulative sum of an array. Like an integral.
	
	dx = dx || 1.0
	
	array  = array.reduce(function(r, a) {
		if (r.length > 0) a += r[r.length - 1];
		r.push(a);
        return r;
		}, []);
	
	return math.multiply(array, dx)
}

function gaussian(x, x0, sigma) {
	// return a normalized gaussian	with the appropriate x0 and sigma
	return math.multiply(math.pow(2*math.PI*math.pow(sigma,2.0), -0.5),math.exp(math.divide(math.subtract(0, math.dotPow(math.subtract(x, x0),2.0)), 2.0*math.pow(sigma, 2.0))));
}

function erf(x, x0, sigma, dx) {
	dx = dx || 1.0
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
	return yf
}


function barrier(V, seq, N2C, laa, Eb, xb, wb, fu) {

	var protseq = seq2charge(seq, N2C);
	var x = math.multiply(math.range(0, protseq.length).toArray(), laa);
	Ue = math.multiply(cumsum(protseq, laa), 1/laa, e*V/kT);
	
	return [x, math.add(Ue, US(x, laa), math.multiply(Eb, erf(x, xb, wb, laa)), math.multiply(fu, x))]
}

function calcEscapeTime2(x, x0, U, D, laa) {
	// calculate passage times
	U = math.subtract(U, U[0]);
	Ur = U.slice().reverse(); // slice is necessary not to reverse in place

	var PiR = math.divide(cumsum(math.exp(U), laa), (math.sum(math.exp(U)) * laa));
	var thetaRiterm = math.multiply(math.dotMultiply(math.exp(U), cumsum(math.dotMultiply(PiR, math.exp(math.multiply(U, -1))), laa)), 1/D);
	var thetaRterm = math.multiply((math.sum(thetaRiterm)*laa),cumsum(math.exp(U), laa), 1/(math.sum(math.exp(U))*laa));
	var thetaR = math.subtract(thetaRterm,cumsum(thetaRiterm, laa));

	var PiL = math.divide(cumsum(math.exp(Ur), laa), (math.sum(math.exp(Ur))*laa)).reverse();
	var thetaLiterm = math.multiply(1/D, math.dotMultiply(math.exp(Ur), cumsum(math.dotMultiply(PiL.slice().reverse(), math.exp(math.multiply(Ur, -1))))));
	var thetaLterm = math.multiply((math.sum(thetaLiterm)*laa), cumsum(math.exp(Ur), laa), 1/(math.sum(math.exp(Ur))*laa));
	var thetaL = math.subtract(thetaLterm, cumsum(thetaLiterm, laa)).reverse();
	
	var pisideL = interp(x0, x, PiL);
	var pisideR = interp(x0, x, PiR);
	
	var tescL = interp(x0, x, math.dotMultiply(thetaL, math.dotPow(PiL, -1)));
	var tescR = interp(x0, x, math.dotMultiply(thetaR, math.dotPow(PiR, -1)));
	
	var tesc = interp(x0, x, math.add(thetaL, thetaR));
	
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
	
	return xsub[ysub.indexOf(math.min(ysub))]
}

var Vs = math.range(25., 80., 5.)
var fu = -0.1
var D = math.pow(10., 5.6)

var tesc = [];
for (var i=0; i<Vs.length; i++) {
	var V = Vs[i];
	var Ub = barrier(V, asynseq, N2C, laa, Eb, xb, wb, fu)
	var x = Ub[0]
	var U = Ub[1]
	x0 = minInRange(x, U, 5., 20.)
	tesc.push(calcEscapeTime2(x, x0, U, D, laa))
}


var protseq = seq2charge(asynseq, false);
var x = math.multiply(math.range(0, protseq.length).toArray(), laa);
var V = 40;

Ub = barrier(V, asynseq, false, laa, 16., 14., 3.5, fu)
document.write(Ub[1])


var Ue = math.multiply(cumsum(protseq, laa), 1/laa, e*V/kT)

//construct potential
var U = math.add(Ue, US(x, laa), math.multiply(16., erf(x, 14., 3.5, laa)), math.multiply(fu, x))

// calculate passage times
U = math.subtract(U, U[0])
Ur = U.slice().reverse() // slice is necessary not to reverse in place

//for (var i=0; i<x.length; i++) {document.write(i, '\t', x[i], '\t', U[i], '<br>');}

var PiR = math.divide(cumsum(math.exp(U), laa), (math.sum(math.exp(U)) * laa));
var thetaRiterm = math.multiply(math.dotMultiply(math.exp(U), cumsum(math.dotMultiply(PiR, math.exp(math.multiply(U, -1))), laa)), 1/D)
var thetaRterm = math.multiply((math.sum(thetaRiterm)*laa),cumsum(math.exp(U), laa), 1/(math.sum(math.exp(U))*laa))
var thetaR = math.subtract(thetaRterm,cumsum(thetaRiterm, laa))

var PiL = math.divide(cumsum(math.exp(Ur), laa), (math.sum(math.exp(Ur))*laa)).reverse()
var thetaLiterm = math.multiply(1/D, math.dotMultiply(math.exp(Ur), cumsum(math.dotMultiply(PiL.slice().reverse(), math.exp(math.multiply(Ur, -1))))))
var thetaLterm = math.multiply((math.sum(thetaLiterm)*laa), cumsum(math.exp(Ur), laa), 1/(math.sum(math.exp(Ur))*laa))
var thetaL = math.subtract(thetaLterm, cumsum(thetaLiterm, laa)).reverse()

    //PiL = -(cumtrapz(numerix.exp(Uf), xf, initial=0.) / trapz(numerix.exp(Uf), x))[::-1]
    //thetaL = (-cumtrapz(numerix.exp(Uf)*cumtrapz(numerix.exp(-Uf)*PiL[::-1]/D, xf, initial=0.), xf, initial=0.) + \
     //           trapz(numerix.exp(Uf)*cumtrapz(numerix.exp(-Uf)*PiL[::-1]/D, xf, initial=0.), xf) * cumtrapz(numerix.exp(Uf), xf, initial=0.)/trapz(numerix.exp(Uf), xf))[::-1]




//for (var i=0; i<x.length; i++) {document.write(i, '\t', x[i], '\t', math.dotMultiply(thetaL, math.dotPow(PiL, -1))[i], '<br>');}
//for (var i=0; i<x.length; i++) {document.write(i, '\t', x[i], '\t', PiR[i], '\t',x[i], '\t', interp(x[i], x, PiR), '\t', x[i+1], '\t', PiR[i+1],  '<br>');}










