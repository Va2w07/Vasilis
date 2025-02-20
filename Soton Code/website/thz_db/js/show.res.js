// This script is for displaying the extraction data and downloading
$(function ()
{
	var saveBut = $("#save");	    		
	var downloadBut = $("#download");
	var timePlotRes = $("#timePlotRes");
	var nPlot = $("#nPlot");
	var dbForm = $('#dataDialogForm');
	var author = 'author';
	var description = 'description';
	var ymin = 0;
	var ymax = 10;
	
	// These options are for the ajax form that sends and recives data from the backend.
    var options = { 
            success:	showResponse,
            error:   	ajaxError,
            type:		'post',       
            dataType:	'json',
            async:		false
    }; 

    var thickSpin = $('#thick').spinner({step:0.01});
    
    // This dialog loads when the tab is loaded.
    // it collects the thickness of the sample in mm from the user and then parses the thz scan data to the backend
	$('#thickDialogForm').dialog({
		autoOpen: true,
		height: 'auto',
		width: 'auto',
		modal: true,
		buttons: {
			"Do extraction": function()
			{
				thickness = thickSpin.spinner("value");
				$('#upthickness').val(thickness);
				$("#upxRVal").val(xRVal);
				$("#upyRVal").val(yRVal);
				$("#upxSVal").val(xSVal);
				$("#upySVal").val(ySVal);
				$("#upData").ajaxForm(options);
				$("#upData").submit();	
				$(this).dialog("close");
			}
		},
		open: function(){
		    $("#thickDialogForm").keypress(function(e) {
		        if (e.keyCode == $.ui.keyCode.ENTER) {
		        	e.preventDefault();
		        }
		    });
		},
		close: function(){}
	});
	
	var timeoptions = 
	{
		series: { lines: { show: true }, shadowSize: 0 },
    	zoom: { interactive: true },
    	pan: { interactive: true }
    };
	
	// Download button for dowloading the complex refractive index as a csv file.
	downloadBut
		.button()
		.click(function( event ){
    		var fname = "refractive_index.csv";
    		var csvContent = "Frequency (THz), n.Real, n.Imag\n";
				
			nReal.forEach(function(infoArray, index){
				dataString = infoArray.join(",") + "," + String(nImag[index][1]);
				csvContent += dataString + "\n";
			});
		
			var pom = document.createElement('a');
			var blob = new Blob([csvContent],{type: 'text/csv;charser=utf-8;'});
			var url = URL.createObjectURL(blob);
			pom.href = url;
			pom.setAttribute('download',fname);
			pom.click();
		});
		
	// if the backend responds extract the data and plot it	 
	function showResponse(data)  { 
		if (data == null){
			return;
		}
		nI = data.nI;
		nI = nI.split(',').map(function(x){return parseFloat(x)});

		nR = data.nR;
		nR = nR.split(',').map(function(x){return parseFloat(x)});

		f = data.f;
		f = f.split(',').map(function(x){return parseFloat(x)});
		
		end = parseFloat(data.end);
		start = parseFloat(data.start);

		doPlot();
	} 

	// error handler if there is a problem from the backend, this needs to be finished with helpfull respnses
	function ajaxError(){
	}
	
	// plots the data
	function doPlot(){	
		df = f[1] - f[0];
		
		nReal = getData(f, nR, 0, f.length);
		nImag = getData(f, nI, 0, f.length);

		ymin = Math.min.apply(null, nI.slice(start,end)) - 1;
		ymax = Math.max.apply(null, nR.slice(start,end)) + 1;

		var noptions = {
				
			series: { lines: { show: true }, shadowSize: 0 },

			xaxis: {min: start*df, max: end*df},
		    yaxis: {min: ymin, max: ymax},
		    
			zoom: { interactive: true },
	    	pan: { interactive: true }

		};

		tP = $.plot(timePlotRes, [dataR, dataS], timeoptions);
		nP = $.plot(nPlot, [nReal, nImag], noptions);

	}

	//get data over the desired range.	
	function getData(x, y, start, finish){
		var d = [];
		for (var i = start; i < finish; i++){
			d.push([x[i], y[i]]);
		}
		return d;
	}
});