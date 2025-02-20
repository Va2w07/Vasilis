<!DOCTYPE html>
<html lang="en">

    <head>
    	<script language="javascript" type="text/javascript" src="js/show.res.js"></script>	
   
   		<form id = "upData" action="/getRefractive" method="post" enctype="multipart/form-data">  
			<input id = "upxRVal" type ="hidden" name = "xRVal" value = "" />
			<input id = "upyRVal" type ="hidden" name = "yRVal" value = "" />
			<input id = "upxSVal" type ="hidden" name = "xSVal" value = "" />		
			<input id = "upySVal" type ="hidden" name = "ySVal" value = "" />
			<input id = "upthickness" type ="hidden" name = "thick" value = "" />
		</form>  
    </head>
    		
    <body>
      	<div id = "timePlotRes" style="height:300px"></div>
      	<div id = "nPlot" style="height:300px"></div>
	    	 
	    <input id = "download" type="submit" name="submit" value="Download Refractive Index" />
	    
	    <form id = "dbformData" action="/addToDB" method="post" enctype="multipart/form-data">  
			<input id = "dbtitle" type = "hidden" name = "title" value = "" />
			<input id = "dbauthor" type = "hidden" name = "author" value = "" />
			<input id = "dbdescription" type = "hidden" name = "description" value = "" />
			<input id = "dbxRVal" type ="hidden" name = "xRVal" value = "" />
			<input id = "dbyRVal" type ="hidden" name = "yRVal" value = "" />
			<input id = "dbxSVal" type ="hidden" name = "xSVal" value = "" />
			<input id = "dbySVal" type ="hidden" name = "ySVal" value = "" />
			<input id = "dbthick" type ="hidden" name = "thick" value = "" />
			<input id = "dbnReal" type ="hidden" name = "nReal" value = "" />
			<input id = "dbnImag" type ="hidden" name = "nImag" value = "" />
			<input id = "dbfreq" type ="hidden" name = "freq" value = "" />				
		</form>  
	    		
		<div id = "thickDialogForm" title = "Thickness">
			<p class = "validateTips"> Please give the thickness of the sample in mm</p>
			<form>
				<fieldset>
	        		<label for = 'thickness'>Thickness (mm)</label>
       				<input type="spinner" name="Thickness" id = 'thick' value = "0.01">	
				</fieldset>
			</form>
		</div>            	          
    </body>
</html>    