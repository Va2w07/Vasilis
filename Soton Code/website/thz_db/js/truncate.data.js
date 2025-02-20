/*
This script is for the truncation screen
*/

$(function ()
{
	var w =0;
	var NN = 0;

	var tstart = xRVal[0];
	var tend = xRVal[xRVal.length-1];

	var dt = Math.abs((xRVal[xRVal.length-1]-xRVal[0])/xRVal.length);
	var tOffset = xSVal[ySVal.indexOf(Math.max.apply(null, ySVal))] - xRVal[yRVal.indexOf(Math.max.apply(null, yRVal))];

	var wnew;
	var W1;
	var W2;

	var rslider = $("#slider-range-ref");
    var sslider = $("#slider-s-ref");

    var rPlot = $("#timeRef");
    var sPlot = $("#timeSam");
	var tabsBar = $("#tabs");

    var s = 100;
    			
    var offset;
    
    //Plotting options for refrence and sample
    var refoptions = {
    	series: { lines: { show: true }, shadowSize: 0 },
    	xaxis: { zoomRange: [0.1, 40], panRange: [0, 40] },
    	yaxis: { zoomRange: [0.1, 40], panRange: [-5, 35] },
    	selection: {
    		mode: "x",
		   	navigate: true
	    }
    };

    var samoptions = {
    	series: { lines: { show: true }, shadowSize: 0 },
    	xaxis: { zoomRange: [0.1, 40], panRange: [0, 40] },
    	yaxis: { zoomRange: [0.1, 40], panRange: [-5, 35] },
    	selection: {
    		mode: "x",
	    	navigate: true
    	}
   	};
	
	// Slider that chooses the range the window has the value of 1
   	rslider.slider({
   		range: true,
	    min: xRVal[0],
	    max: xRVal[xRVal.length-1],
	    values: [tstart/dt, tend/dt],
	    slide: function( event, ui) 
	    {
	    	tstart = ui.values[0];
	    	tend = ui.values[1];
	    	offset = Math.round(xRVal.length - tstart/dt - (tend-tstart)/(2*dt)); 
	    	
	    	//trucWindow creates the truncation wind with the corrent shape we want	however it needs to be shifted
	    	// for both the sample and refrence	    				
	    	wnew = truncWindow(s, Math.round((tend-tstart)/dt));

	    	// Shift the window relative to the position of the peak in time for both data.
    		W1 = shiftW(wnew, offset, dataR);
    		W2 = shiftW(wnew, Math.round(offset-tOffset/dt), dataS);
	    	
	    	//Set window data for plot	
    		rP.setData([dataR, W1]);
    		sP.setData([dataS, W2]);
			
			//Draw window
    		rP.draw();
    		sP.draw();
	    	
	    	// Draw the higlighted selected area.
    		rP.setSelection({xaxis:{from: tstart, to: tend}});
    		sP.setSelection({xaxis:{from: tstart+tOffset, to: tend+tOffset}});
		},
		step: dt
	});

	// This slide changes the slope of the cosine of the truncation windoww    		
   	sslider.slider({
  	   	range: "min",
   	   	min: 0,
   	   	max: (xRVal[xRVal.length-1]-xRVal[0])/2, // Max is 50% of the data
   	   	values: [s*dt],
   	   	slide: function( event, ui) {
   	   		s = ui.values[0]/dt;
   	   		//Build the truncation window
   			wnew = truncWindow(s, Math.round((tend-tstart)/dt));

   			//Shift the truncation window for refrence and sample to it is in the same place relative to the peak
   			W1 = shiftW(wnew, offset, dataR);
	   		W2 = shiftW(wnew, Math.round(offset-tOffset/dt), dataS);
	    	
	    	//Set window data for plot
	   	    rP.setData([dataR, W1]);
	   	    sP.setData([dataS, W2]);
	    	
	    	//Draw the window. 
	   	    rP.draw();
	   	    sP.draw();
	   	    // As the width of the main box did not change we do not need to re-render the highlighted box.
   	   	},
	    step: dt
	});
		
   	// This function builds the truncation window which is a taoered cosine (Turky Window)
	function truncWindow(s,  NN){
		// Set the window to be double that of the data so it can be rotated without looping back 
		// and being applied to data we do not want it to touch
		var w = Array.apply(null, new Array(xRVal.length*2)).map(Number.prototype.valueOf,0); 

		// sets the uper and lower bounds of where w = 1
		var lower = Math.round(w.length/2 - NN/2); 
		var upper = Math.round(w.length/2 + NN/2);
		
		// this builds the truncation window
		for (var i = Math.round(w.length/2 - NN/2 - s); i < Math.round(w.length/2 + NN/2 + s); i++) {
			if (i < lower) {
				w[i] = (Math.cos(Math.PI*((i+NN/2+s-w.length/2)/s + 1)) + 1)/2; 
			}
			else if (i < upper) {
				w[i] = 1;
			}
			else {
				w[i] = (Math.cos(Math.PI*(i-NN/2-w.length/2)/s ) + 1)/2;
			}
		}				
		return w;
	}
	
	// This function shifts the window to align it with the thz peaks
	function shiftW(w,  offset,  data) {
		var wArray = [];
		offt = Math.round(data[0][0]/dt);
		for (var i = 0; i < w.length/2; i++) {
			wArray.push([data[i][0], w[(i+offset+offt) % w.length]]);					
		}
		
		return wArray
	}

	// This function performes the truncation to the data to be parsed
	function truncData(){
    	for (var i = 0; i<xSVal.length; i++){
			dataR[i][1] = yRVal[i] * W1[i][1]; 
			dataS[i][1] = ySVal[i] * W2[i][1]; 
    	}
	}

	// Next button for going to the extraction
	$("#truncBut")
		.button()
		.click(function () {
			truncData(); // Truncate the data to be parsed
			tabsBar
				.tabs(
				{
					disabled: [1], // Disable truncation tab
					active: [2] // Activate the Extraction tab
				});
		});

	// Plot the data
   	rP = $.plot(rPlot, [dataR], refoptions);
   	sP = $.plot(sPlot, [dataS], samoptions);
    
    // set up slider
   	rslider.css('left', rP.getPlotOffset().left);    		
   	rslider.width(rslider.width() - (rP.getPlotOffset().right + rP.getPlotOffset().left));
    
    // Set selection area
   	rP.setSelection({xaxis: {from: xRVal[0], to: xRVal[xRVal.length -1]}});
	sP.setSelection({xaxis:{from: xSVal[0], to: xSVal[xSVal.length -1]}});
});