<!DOCTYPE html>
<html lang="en">
    <head>
        <meta charset="utf-8">
    </head>

    <script language="javascript" type="text/javascript" src="js/truncate.data.js"></script>	

    <body>
	    <div id = "timeRef" style=" height:300px"></div>
		<div id="slider-range-ref" class = "ui-slider ui-slider-horizontal ui-widget ui-widget-content ui-corner-all"></div>	
		<div id="slider-s-ref" class = "ui-slider ui-slider-horizontal ui-widget ui-widget-content ui-corner-all"></div>	
	    <div id = "timeSam" style="height:300px"></div>
		<div id="slider-range-sam"></div>	
	    <input id = "truncBut" type="submit" name="submit" value="Next" />	    	    
	</body>
	
	<form id = "formData" action="/doExtraction" method="post" enctype="multipart/form-data">  
		<input id = "title" type = "hidden" name = "title" value = "" />
		<input id = "thickness" type ="hidden" name = "thickness" value = "" />
		<input id = "xRVal" type ="hidden" name = "xRVal" value = "" />
		<input id = "yRVal" type ="hidden" name = "yRVal" value = "" />
		<input id = "xSVal" type ="hidden" name = "xSVal" value = "" />
		<input id = "ySVal" type ="hidden" name = "ySVal" value = "" />
		<input id = "offset" type ="hidden" name = "offset" value = "" />
		<input id = "width" type ="hidden" name = "width" value = "" />
		<input id = "w" type ="hidden" name = "w" value = "" />
	</form>  
</html>