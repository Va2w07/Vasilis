<!DOCTYPE html>
<html lang="en">
    <head>
        <meta charset="utf-8">		
	</head>
	<script language="javascript" type="text/javascript" src="js/upload.data.js"></script>	

	<body>
		<div class="ui-widget">	
    		<label for = 'reference'>Reference Data</label>
        	<input id = 'reference' class="button ui-widget-content ui-corner-all" type="file" name="file[]"/>

    		<label for = 'sample'>Sample Data</label>
        	<input id = 'sample' class="button ui-widget-content ui-corner-all" type="file" name="file[]"/>
        	
	        <label for = "time">Time Unit Of Your Data</label>       	
    		<select name="Time" id = "time">
        		<option value = '1'>s</option>
        		<option value = '1e-12'>ps</option>
        	</select>        
        	<div id = "timePlot" style="height:500px"></div>
	    	<p class="timemessage"></p>
	    	
	    	<button type = "submit" id = "clearPlot" name = "clearPlot">Clear Data</button>	    	
	    	<button type = "submit" id = "next" name = "next">Next</button>	    	
        	   		
		</div>

	</body>
	
	<div id = "getSepForm" title = "Please select the seporator used in you file ">
		<p class = "validateTips"> Data must be in two columns where the first column is time and the second the time scan.</p>
		<form>
	        <label for = "seporator">Separator</label>       	
	 		<select name="seporator" id = "seporator">
       			<option value = '	'>Tab</option>
       			<option value = ','>Comma</option>
       			<option value = ' '>Space</option>
       		</select>
       		<input name = "refOrsam" type = "hidden" id = "refOrSam">        					
		</form>
	</div>          
	
	<div id = "doTrunDialog" title = "Would you like to truncate the data">
		<p class = "validateTips"> If your sample scan has reflections then truncating the data to remove them will reduce oscillations in the refractive index...</p>
	</div>          
</html>