/*
This script is for uploading and displaying the raw thz tds scans from the user.
*/

$(function()
{

	var timeCombo = $('#time');
	var refUp = $('#reference');
	var samUp = $('#sample');
	var timePlot = $("#timePlot");
	var sepCombo = $("#seporator");
	var sepForm = $("#getSepForm");
	var refOrsam = $("#refOrSam");
	var clrPlot = $("#clearPlot");
	var trunDialog = $("#doTrunDialog");
	var nextBut = $('#next');
	var tabsBar = $("#tabs");
	
	var refFile;
	var samFile;
	var thickness;

	// Reset values to their defaults
	resetVals();    
	
	nextBut
		.button()
		.click(function(event)
		{
			// copy data to refrence and sample arrays.
			for (var i = 0; i < dataR.length; i++) 
			{
				xRVal[i] = dataR[i][0];
				yRVal[i] = dataR[i][1];
			}

			for (var i = 0; i < dataS.length; i++)
			{
				xSVal[i] = dataS[i][0];
				ySVal[i] = dataS[i][1];
			}

			// Error Catchers for common problems.
			if(xRVal.length == 0 || yRVal.length == 0){
				alert('Data cannot be empty.');
				return;
			} 

			if(xRVal.length != yRVal.length){
				alert('Time data must have the same number of points.')
				return;
			} 

			if( Math.abs((xRVal[1] - xRVal[0]) - (xSVal[1] - xSVal[0])) > 0.1*(xSVal[1] - xSVal[0])){
				alert('Time step must be the same for both scans')
				return;				
			} 
			
			// open dialog to ask if we want to do truncation and go to the next step
			trunDialog.dialog("open");
		});
	
	// Truncation dialog to load when next is clicked
	trunDialog
		.dialog({
			autoOpen: false,
			height: 'auto',
			width: 'auto',
			modal: true,
			buttons: 
			{
				Yes: function()
				{
					tabsBar
						.tabs(
						{
							disabled: [2], //Disable extraction tab if open
							active: [1] //Activate Truncation tab
						});
					$(this).dialog('close');
				},
				No: function()
				{
					tabsBar
						.tabs(
						{
							disabled: [1], // Disable truncation tab if open
							active: [2] // Activate extraction tab
						});

					$(this).dialog('close');
				}
			},
		});
	
	// Clears the plots
	clrPlot
		.button()
		.click(function(event){
			// We want to reset to defaults and null the uploaded files
			resetVals();
			refUp.val(null);
			samUp.val(null);
			
			tP = $.plot(timePlot, [[]], timeoptions);
			
			tabsBar
				.tabs(
				{
					disabled: [1,2] // Disable truncation and extraction tabs
				});
	    });
	
	//This dialog gets the datas delimiter, comma, tab or space and then calls the extracion.
	sepForm.dialog({
		autoOpen: false,
		height: 'auto',
		width: 'auto',
		modal: true,
		buttons: {
			"Upload Data": function(){
				switch (parseInt(refOrsam.val())) // 1 for ref scan, 2 for sample scan
				{
				case 1:
					extractData( refFile, 1, sepCombo.val() );
					break;
				case 2:
					extractData( samFile, 2, sepCombo.val() );
					break;
				}
				
				$(this).dialog("close");
			}
		},
		Cancel: function(){
			$(this).dialog("close");
		},
		close: function(){
		}
	});

	sepCombo
		.combobox();
	
	// This combo box is for selecting the time step unit (s or ps)
	timeCombo
		.val(tUnit) //Set unit to current time unit or defualt.
		.combobox(
		{
			select: function(request, response)
				{
					tUnit = timeCombo.val();
					for (var i = 0; i < dataR.length; i++)
					{
						dataR[i][0] = timeOldR[i] * tUnit; 
					}

					for (var i = 0; i < dataS.length; i++)
					{
						dataS[i][0] = timeOldS[i] * tUnit; 
					}
					tP = $.plot(timePlot, [dataR, dataS], timeoptions);
				}
		});

	//Upload buttons for refrence and sample.

	refUp
		.bind('change', function (evt)
		{
			if (evt.target.files.length == 0)
			{
				evt.target.files = refFile;
			}
			else
			{
				refOrsam.val(1);
				refFile = evt.target.files;
				sepForm.dialog("open");
			}
		});

	samUp
		.bind('change', function (evt)
		{
			if (evt.target.files.length == 0)
			{
				evt.target.files = samFile;
			}
			else
			{
				refOrsam.val(2);
				samFile = evt.target.files;
				sepForm.dialog("open");
			}
		});

	// plot options
	var timeoptions = 
	{
		series: { lines: { show: true }, shadowSize: 0 },
    	zoom: { interactive: true },
    	pan: { interactive: true }
    };

	var normflag = false;	
	var normalise = 1;
	var snorm = 1;
	var offTime = 0;
	
	//This function extracts the data from the upload files, it is called by the upload buttons
	// This funcation is a bit of a mess and can definatly be trimmed, sorry.
    function extractData(file, dataSet, sep)
    {
    	var reader = new FileReader();
    	reader.readAsText(file[0]);
    	
    	switch(dataSet) // 1 for ref, 2 for sample.
    	{
    		case 1:
    			reader.onload = function(event)
    			{
    				var csv = event.target.result;
    				var newData = $.csv.toArrays(csv, {separator: sep});

    				dataR = newData;
    				timeOldR = [];
    				
    				normalise = getMaxY(dataR);
    				
    				offTime = dataR[0][0] * 1e12 * tUnit;
    				
	   				for (var i = 0; i < newData.length; i++)
	   				{
	   					timeOldR[i] = dataR[i][0] * 1e12; 
	   					dataR[i][0] = timeOldR[i] * tUnit - offTime; 
	   					dataR[i][1] = parseFloat(dataR[i][1])/normalise; 
	   				}
	   				
	   				if (dataS.length != 0)
	   				{
   						for (var i = 0; i < dataS.length; i++)
	   					{
   							dataS[i][0] = dataS[i][0] - offTime;
	   						dataS[i][1] = dataS[i][1] * snorm/normalise;
	   					}
	   				}
	   				tP = $.plot(timePlot, [dataR, dataS], timeoptions);    	
	   				normflag = true;
    			};
    			break;
    	    
    		case 2:	
    			reader.onload = function(event)
    			{
    				var csv = event.target.result;
    				var newData = $.csv.toArrays(csv, {separator: sep});

    				dataS = newData;	
    				timeOldS = [];

    				if (normflag == false)
    				{
        				normalise = getMaxY(dataS);
    				}
    				
    				for (var i = 0; i < dataS.length; i++)
	   				{
	   					timeOldS[i] = dataS[i][0] * 1e12;
	   					dataS[i][0] = timeOldS[i] * tUnit - offTime;
	   					dataS[i][1] = parseFloat(dataS[i][1])/normalise;	   					
	   				}
	   				tP = $.plot(timePlot, [dataR, dataS], timeoptions);
	   				snorm = normalise;
	   				normalise = 1;
	   				normflag = false;
	    		};
				break;
		}
    	reader.onerror = function(){alert('Cannot read file')};
    }

    //Gets the maximum value so we can normalise the data.    
    function getMaxY(dataArray)
    {
    	var maxVal = 0;

    	for (var i = 0; i  < dataArray.length; i++)
    	{
    		if (parseFloat(dataArray[i][1]) > maxVal)
    		{
    			maxVal = parseFloat(dataArray[i][1]);
    		}
    	}
    	return maxVal;
    }

    //resets all the globals 
    function resetVals()
	{
		dataR = [];
		dataS = [];

		timeOldR = [];
		timeOldS = [];

		tUnit = 1;

		nI = [];
		nR = [];
		f = [];
		start = 0;
		end = 0;
		thickness = 0;

		nReal = [];
		nImag = [];

		df = 0;

		xRVal = [];
		yRVal = [];

		xSVal = [];
		ySVal = [];

		refUp.val(null);
		samUp.val(null);
		
		tP = $.plot(timePlot, [[]], timeoptions);
			
		tabsBar
			.tabs(
			{
				disabled: [1,2]
			});
	}

	tP = $.plot(timePlot, [dataR, dataS], timeoptions);
});