/*
This script manages the tabs on the top of the page.
0 - Upload tab
1 - Truncation tab
2 - 
*/
$(function() 
{
	var tUnit = 1; // Defualt time unit is 1
	var disabledTabs = [1,2]; // By defualt only the upload tab (0) is active
	var tabsBar = $("#tabs");

	tabsBar
		.tabs(
		{
			disabled: disabledTabs, 
			beforeLoad: function( event, ui ) 
			{
				ui.jqXHR.error(function() 
				{
					ui.panel.html("Couldn't load this tab. We'll try to fix this as soon as possible.");
				});
			}
		});
});
