/*
	PROTEIN THERMODYNAMICS SIMULATIONS
	Copyright (C) 2021 Johan Pääkkönen, Juha Rouvinen, University of Eastern Finland
	
	This program is free software; you can redistribute it and/or
	modify it under the terms of the GNU General Public License
	as published by the Free Software Foundation; either version 2
	of the License, or (at your option) any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.
	
	You should have received a copy of the GNU General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

/*
	The GNU General Public Licence, version 2, is available at:
	(1) the file LICENSE in this folder
	(2) https://www.gnu.org/licenses/old-licenses/gpl-2.0.html
*/

"use strict";

var S_0 = 1;
var E_0 = 1;
var Q_0 = 1;
var K_D = 1;
var K_D2 = 1;

var scale_absolute = 1;
var xscale_alternative = false;

var fwidth_old = 0;
var fheight_old = 0;

var xaxistype = 0;
var yaxistype = 0;
var axistype_lin = 0;
var axistype_log = 1;
var axistype_time = 2;

var curves = [];
var labels = [];
var legends = [];
var colours = [];
var pie = [];
var datapoints = [];

var xmin = 0, xmax = 0;
var ymin = 0, ymax = 0;

var appmode = -1;
var appmode_listonly  = -1;
var appmode_ligands   = 0;
var appmode_receptors = 1;
var appmode_ligand    = 2;
var appmode_homodimer = 3;

// Set to true to display indicators of the pie segments.
// I made it only for preparing figures for my PhD thesis,
// so do not complain if it behaves weirdly.
var extpiemode = false;

var svg_xmlns = "http://www.w3.org/2000/svg";
var xml_xmlns = "http://www.w3.org/XML/1998/namespace";

var pages = [];
pages[appmode_ligands]   = {title: "Competing ligands simulation",        link: "ligands.htm"   };
pages[appmode_receptors] = {title: "Competing receptors simulation",      link: "receptors.htm" };
pages[appmode_ligand]    = {title: "Ligand binding simulation",           link: "ligand.htm"    };
pages[appmode_homodimer] = {title: "Protein homodimerisation simulation", link: "homodimer.htm" };

function preinit(mode)
{
	appmode = mode;
	
	var linkarea = document.getElementById("linkarea");
	var credits = document.getElementById("credits");
	var brnode, anode, textnode;
	var first = true;
	var i;
	
	credits.innerHTML = "Johan P&auml;&auml;kk&ouml;nen, Juha Rouvinen<br />University of Eastern Finland, 2021";
	
//	credits.innerHTML = "Johan P&auml;&auml;kk&ouml;nen, Juha Rouvinen<br />University of Eastern Finland, 2021" +
//		"<br/><br />If you use this in your own research, please cite our article: <a href=\"\">https://doi.org/XXX.XXX/XXX</a>";
	
	for(i = 0; i < pages.length; i++)
	{
		if(appmode === i || pages[i] === undefined) continue;
		
		if(!first)
		{
			brnode = document.createElement("br");
			linkarea.appendChild(brnode);
			brnode = document.createElement("br");
			linkarea.appendChild(brnode);
		}
		
		first = false;
		
		textnode = document.createTextNode(pages[i].title);
		anode = document.createElement("a");
		anode.appendChild(textnode);
		anode.href = pages[i].link;
		linkarea.appendChild(anode);
	}
	
	if(!(appmode >= 0)) return;
	
	var datatable = document.getElementById("datatable");
	
	for(i = -1; i < datalabels[appmode].length; i++)
	{
		var ele = document.createElement("tr");
		var td1 = document.createElement(i === -1 ? "th" : "td");
		var td2 = document.createElement(i === -1 ? "th" : "td");
		var td3 = document.createElement(i === -1 ? "th" : "td");
		var td4 = document.createElement(i === -1 ? "th" : "td");
		
		if(i === -1)
		{
			td1.innerHTML = "Species";
			td2.innerHTML = "Conc. (mol l<sup>&minus;1</sup>)";
			td3.innerHTML = "Conc. (g l<sup>&minus;1</sup>)";
			td4.innerHTML = "Proportion (%)";
		}
		else
		{
			td1.innerHTML = datalabels[appmode][i].substring(1, datalabels[appmode][i].length - 1);
			td2.innerHTML = "0";
			td3.innerHTML = "0";
			td4.innerHTML = "0";
			td1.setAttribute("id", "data" + i + "l");
			td2.setAttribute("id", "data" + i + "a");
			td3.setAttribute("id", "data" + i + "b");
			td4.setAttribute("id", "data" + i + "c");
		}
		
		ele.appendChild(td1);
		ele.appendChild(td2);
		ele.appendChild(td3);
		ele.appendChild(td4);
		datatable.appendChild(ele);
	}
	
	var calcinfo = document.getElementById("calcinfo");
	
	calcinfo.innerHTML =
		"Give a list of pairs of values: <strong>x-value\xA0y-value</strong>. The curve of your choice will be least-squares fitted<br />" +
		"to the data by finding the optimal slider position of the dissociation constant.<br />" +
		"The number of data points must be greater than or equal to the number of free parameters.<br />" +
		"Note that this is not a proper analytical tool. Results are imprecise and possibly inaccurate.";
	
	var figure_div = document.getElementById("figure_div");
	
	figure_div.style.width = "625px";
	figure_div.style.minWidth = "450px";
	figure_div.style.height = "450px";
	figure_div.style.minHeight = "300px";
	
	/*
	figure_div.style.width = "585px";
	figure_div.style.minWidth = "585px";
	figure_div.style.height = "300px";
	figure_div.style.minHeight = "300px";
	*/
	
	figure_div.innerHTML = 	
		"<svg xmlns=\"http://www.w3.org/2000/svg\" id=\"figure_svg\" width=\"96%\" height=\"100%\" viewBox=\"0 0 600 450\" preserveAspectRatio=\"none\" style=\"font-family: Arial, Helvetica, sans-serif; font-size: 12px\">\n" +
			"<rect x=\"0\" y=\"0\" width=\"100%\" height=\"100%\" style=\"fill: white; stroke: none;\" />\n" +
			"<g id=\"svg_curvegroup\" style=\"fill: none; stroke-width: 3px; stroke-linejoin: round; stroke-linecap: butt;\"></g>\n" +
			"<g id=\"svg_datapointgroup\"></g>\n" +
			"<g id=\"svg_overlaygroup\"></g>\n" +
			"<g id=\"svg_linegroup\" style=\"fill: none; stroke: black; stroke-width: 3px;\"></g>\n" +
			"<g id=\"svg_legendgroup\"></g>\n" +
			"<g id=\"svg_piegroup\"></g>\n" +
			"<g id=\"svg_scribblegroup\"></g>\n" +
		"</svg>";
	
	var toggle_calcdiv_checkbox = document.getElementById("toggle_calcdiv_checkbox")
	if(toggle_calcdiv_checkbox)
	{
		// crude hack, but it works...
		// because the function will reverse this back immediately
		toggle_calcdiv_checkbox.checked = !toggle_calcdiv_checkbox.checked;
		toggle_calcdiv();
	}
}

function init()
{
	var i;
	
	if(!(appmode >= 0)) return;
	
	var sliders = document.getElementsByClassName("slider");
	
	for(i = 0; i < sliders.length; i++)
	{
		slider_input(Number(sliders[i].id.substr(6)), true);
	}
	
	var radio1 = document.getElementById("radio1");
	scale_absolute = (radio1 && radio1.checked) ? 1 : 0;
	
	var radio5 = document.getElementById("radio5");
	if(radio5)
	{
		var radio6 = document.getElementById("radio6");
		var radio7 = document.getElementById("radio7");
		
		if(radio5.checked) scale_absolute = 1;
		if(radio6.checked) scale_absolute = 0;
		if(radio7.checked) scale_absolute = 2;
	}
	
	var radio3 = document.getElementById("radio3");
	xscale_alternative = !(radio3 && radio3.checked);
	
	radio_input(xscale_alternative ? 4 : 3, true);
	
	var figure_div = document.getElementById("figure_div");
	fwidth_old = extract_px(figure_div.style.width);
	fheight_old = extract_px(figure_div.style.height);

	document.body.addEventListener("mousemove", mouse_moved);
	document.body.addEventListener("mouseup", mouse_released);
	
	var databox = document.getElementById("databox");
	if(databox) databox.addEventListener("blur", data_changed);
	
	var calcbutton = document.getElementById("calcbutton");
	if(calcbutton) calcbutton.addEventListener("click",
		function()
		{
			var label = document.getElementById("calcstatus");
			if(label)
			{
				label.innerHTML = "Calculating, please wait...";
				label.style.color = "blue";
			}
			
			// This is so that the label has time to change
			// before the code enters a heavy calculation.
			window.setTimeout(calculate_lsf, 100);
		}
	);
	
	document.getElementById("export_button").addEventListener("click", export_svg);
	
	data_changed(false);
	update();
}

function render_text_svg(ele, str)
{
	var fontfamily = "";
	var normalfont = "12px";
	var scriptfont = "10px";
	
	var supbl = -4;
	var subpl = 2.5;
	
	var nextchar = 0;
	var prevchar = -1;
	var search = 0;
	var token = "";
	var baseline = 0;
	var c;
	
	delete_all_children(ele);
	
	do
	{
		search = str.substr(prevchar + 1).search(/[\uEEE0-\uEEEF]/);
		nextchar = prevchar + 1 + search;
		token = str.substring(prevchar + 1, (search === -1) ? undefined : nextchar);
		
		c = (prevchar === -1 ? 0xEEE0 : str.charCodeAt(prevchar));
		
		function new_textspan(bold, italic, script, align, content)
		{
			if(!content) return;

			var span = document.createElementNS(svg_xmlns, "tspan");
			if(bold) span.style.fontWeight = "bold";
			if(italic) span.style.fontStyle = "italic";
			if(script) span.style.fontSize = "10px";
			if(align - baseline !== 0) span.setAttribute("dy", align - baseline);
			baseline = align;

			var textnode = document.createTextNode(content);
			span.appendChild(textnode);

			ele.appendChild(span);
		}

		switch(c)
		{
			case 0xEEE0: // normal
			{
				new_textspan(false, false, false, 0, token);
				break;
			}
			case 0xEEE1: // superscript
			{
				new_textspan(false, false, true, supbl, token);
				break;
			}
			case 0xEEE2: // subscript
			{
				new_textspan(false, false, true, subpl, token);
				break;
			}
			case 0xEEE4: // italic
			{
				new_textspan(false, true, false, 0, token);
				break;
			}
			case 0xEEE5: // italic superscript
			{
				new_textspan(false, true, true, supbl, token);
				break;
			}
			case 0xEEE6: // italic subscript
			{
				new_textspan(false, true, true, subpl, token);
				break;
			}
			case 0xEEE8: // bold normal
			{
				new_textspan(true, false, false, 0, token);
				break;
			}
			case 0xEEE9: // bold superscript
			{
				new_textspan(true, false, true, supbl, token);
				break;
			}
			case 0xEEEA: // bold subscript
			{
				new_textspan(true, false, true, subpl, token);
				break;
			}
			case 0xEEEC: // bold italic
			{
				new_textspan(true, true, false, 0, token);
				break;
			}
			case 0xEEED: // bold italic superscript
			{
				new_textspan(true, true, true, supbl, token);
				break;
			}
			case 0xEEEE: // bold italic subscript
			{
				new_textspan(true, true, true, subpl, token);
				break;
			}
			default:
			{
				prevchar = nextchar;
				continue;
			}
		}
		
		prevchar = nextchar;
	}
	while(search !== -1);
}

function toggle_calcdiv()
{
	var ele1 = document.getElementById("toggle_calcdiv_symbol");
	var ele2 = document.getElementById("calcdiv");
	var ele3 = document.getElementById("toggle_calcdiv_checkbox");
	
	if(!ele3.checked)
	{
		ele3.checked = true;
		ele2.style.display = "block";
		ele1.innerHTML = "\u25B2";
	}
	else
	{
		ele3.checked = false;
		ele2.style.display = "none";
		ele1.innerHTML = "\u25BC";
	}
}

function replace_minus_signs(str)
{
	return str.replace(/\u002D/g, "\u2212");
}

function expval(v, e, minexp, maxexp)
{
	return Math.pow(e, minexp + v * (maxexp - minexp) / 240);
}

var expparams = [
	/*  0 */ [0, 0, 0],
	/*  1 */ [10, -1, 5],
	/*  2 */ [10, -7, -1],
	/*  3 */ [10, -9, -1],
	/*  4 */ [10, 0, 6],
	/*  5 */ [10, -9, -1],
	/*  6 */ [10, -5, -1],
	/*  7 */ [10, -9, -1],
	/*  8 */ [0, 0, 0],
	/*  9 */ [10, -9, -1],
	/* 10 */ [10, -9, -1],
	/* 11 */ [10, -9, 1],
	/* 12 */ [10, -5, 5],
	/* 13 */ [10, -11, -1],
	/* 14 */ [10, -12, -4],
	/* 15 */ [10, -12, -4]
];

function expval_wrap(v, index)
{
	return expval(v, expparams[index][0], expparams[index][1], expparams[index][2]);
}

function slider_input(index, noupdate)
{
	var slider = document.getElementById("slider" + index);
	var label = document.getElementById("value" + index);
	
	var rawval = (!slider ? undefined : Number(slider.value)); // [0,240]
	var val = 0;
	var valstr = "";
	var unitstr = "";
	
	switch(index)
	{
		case 3: // [S]_0
		{
			val = S_0 = expval_wrap(rawval, index);
			unitstr = "mol\xA0l\uEEE1\u22121\uEEE0";
			break;
		}
		case 5: // [E]
		{
			val = E_0 = expval_wrap(rawval, index);
			unitstr = "mol\xA0l\uEEE1\u22121\uEEE0";
			break;
		}
		case 7: // K_D
		{
			val = K_D = expval_wrap(rawval, index);
			unitstr = "mol\xA0l\uEEE1\u22121\uEEE0";
			break;
		}
		case 8: // delta_G
		{
			val = -8.31446261815324 * 298.15 * Math.log(K_D) * 1e-3;
			unitstr = "kJ\xA0mol\uEEE1\u22121\uEEE0";
			break;
		}
		case 9: // K_D 2
		{
			val = K_D2 = expval_wrap(rawval, index);
			unitstr = "mol\xA0l\uEEE1\u22121\uEEE0";
			break;
		}
		case 10: // [E]
		{
			val = Q_0 = expval_wrap(rawval, index);
			unitstr = "mol\xA0l\uEEE1\u22121\uEEE0";
			break;
		}
		case 16: // K_A (association constant)
		{
			val = 1.0 / K_D;
			unitstr = "l\xA0mol\uEEE1\u22121\uEEE0";
			break;
		}
		case 17: // delta_G 2
		{
			val = -8.31446261815324 * 298.15 * Math.log(K_D2) * 1e-3;
			unitstr = "kJ\xA0mol\uEEE1\u22121\uEEE0";
			break;
		}
		case 18: // K_A 2 (association constant)
		{
			val = 1.0 / K_D2;
			unitstr = "l\xA0mol\uEEE1\u22121\uEEE0";
			break;
		}
	}
	
	var labelstr = "";
	
	var abs = Math.abs(val);
	var magnitude = Math.floor(Math.log10(abs));
	var sign = (val < 0) ? "\u2212" : "";
	
	if(magnitude === 0 || index === 17 || index === 8)
	{
		valstr = sign + abs.toFixed(2);
	}
	else
	{
		valstr = sign + (abs / Math.pow(10, magnitude)).toFixed(2) + "\xA0\u22C5\xA0" + "10\uEEE1" + magnitude + "\uEEE0";
	}
		
	labelstr = replace_minus_signs(valstr) + "\xA0" + unitstr;
	
	while(label.hasChildNodes()) label.removeChild(label.firstChild); 
	
	var nextchar = 0;
	var prevchar = -1;
	var search = 0;
	var elementnode = null;
	var textnode = null;
	
	do
	{
		search = labelstr.substr(prevchar + 1).search(/[\uEEE0\uEEE1\uEEE2]/);
		nextchar = prevchar + 1 + search;
		textnode = document.createTextNode(labelstr.substring(prevchar + 1, (search === -1) ? undefined : nextchar));
		
		if(prevchar === -1 || labelstr.charCodeAt(prevchar) === 0xEEE0)
		{
			label.appendChild(textnode);
		}
		else if(labelstr.charCodeAt(prevchar) === 0xEEE1)
		{
			elementnode = document.createElement("sup");
			elementnode.appendChild(textnode);
			label.appendChild(elementnode);
		}
		else if(labelstr.charCodeAt(prevchar) === 0xEEE2)
		{
			elementnode = document.createElement("sub");
			elementnode.appendChild(textnode);
			label.appendChild(elementnode);
		}
		
		prevchar = nextchar;
	}
	while(search !== -1);
	
	if(index === 7)
	{
		slider_input(8, true);
		slider_input(16, true);
	}
	
	if(index === 9)
	{
		slider_input(17, true);
		slider_input(18, true);
	}
	
	if(!noupdate) update();
}

function radio_input(index, noupdate)
{
	switch(index)
	{
		case 1: // absolute y-axis
		{
			scale_absolute = 1;
			break;
		}
		case 2: // relative y-axis
		{
			scale_absolute = 0;
			break;
		}
		case 3: // x-axis type 1
		{
			xscale_alternative = 0;
			
			if(appmode === appmode_ligand)
			{
				document.getElementById("title3").title = "Total ligand concentration";
				document.getElementById("label3_init").style.display = "inline-block";
				document.getElementById("label3_free").style.display = "none";
			}
			else if(appmode === appmode_homodimer)
			{
				document.getElementById("title5").title = "Total protein concentration";
				document.getElementById("label5_init").style.display = "inline-block";
				document.getElementById("label5_free").style.display = "none";
			}
			else if(appmode === appmode_ligands)
			{
				document.getElementById("fixval5").disabled = true;
				document.getElementById("fixval10").disabled = false;
			}
			
			if(appmode === appmode_receptors)
			{
				document.getElementById("title5").title = "Total ligand concentration";
				document.getElementById("label5_init").style.display = "inline-block";
				document.getElementById("label5_free").style.display = "none";
			}
			
			break;
		}
		case 4: // x-axis type 2
		{
			xscale_alternative = 1;
			
			if(appmode === appmode_ligand)
			{
				document.getElementById("title3").title = "Free ligand concentration";
				document.getElementById("label3_init").style.display = "none";
				document.getElementById("label3_free").style.display = "inline-block";
			}
			else if(appmode === appmode_homodimer)
			{
				document.getElementById("title5").title = "Free protein concentration";
				document.getElementById("label5_init").style.display = "none";
				document.getElementById("label5_free").style.display = "inline-block";
			}
			else if(appmode === appmode_ligands)
			{
				document.getElementById("fixval5").disabled = false;
				document.getElementById("fixval10").disabled = true;
			}
			
			if(appmode === appmode_receptors)
			{
				document.getElementById("title5").title = "Free ligand concentration";
				document.getElementById("label5_init").style.display = "none";
				document.getElementById("label5_free").style.display = "inline-block";
			}
			
			break;
		}
		case 11: // x-axis type 3
		{
			xscale_alternative = 2;
			
			if(appmode === appmode_receptors)
			{
				document.getElementById("fixval5").disabled = false;
				document.getElementById("fixval10").disabled = true;
			}
			
			break;
		}
		case 5: // absolute logarithmic y-axis
		{
			scale_absolute = 1;
			break;
		}
		case 6: // absolute linear y-axis
		{
			scale_absolute = 0;
			break;
		}
		case 7: // specificity axis
		{
			scale_absolute = 2;
			break;
		}
		default:
		{
			break;
		}
	}
	
	if(!noupdate) update();
}

function remove_from_parent(ele)
{
	if(!ele) return null;
	return ele.parentNode.removeChild(ele);
}

function delete_all_children(ele)
{
	if(!ele) return;
	while(ele.hasChildNodes()) ele.removeChild(ele.firstChild);
}

var resizing = false;

function mouse_moved()
{
	var figure_div = document.getElementById("figure_div");
	var fwidth = extract_px(figure_div.style.width);
	var fheight = extract_px(figure_div.style.height);
	
	if(!(fwidth === fwidth_old && fheight === fheight_old))
	{
		fwidth_old = fwidth;
		fheight_old = fheight;
		resizing = true;
		update(false);
	}
}

function mouse_released()
{
	if(resizing)
	{
		resizing = false;
		update();
	}
}

function data_changed(do_update)
{
	var databox = document.getElementById("databox");
	var lines = [databox.value];
	var n, p, num1, num2;
	var line, part1, part2;
	var datastr = "";
	
	datapoints = [];
	
	while(1)
	{
		n = lines.length - 1;
		p = lines[n].indexOf("\n");
		
		if(p === -1) break;
		
		part1 = lines[n].substring(0, p);
		part2 = lines[n].substring(p + 1);
		
		lines[n] = part1;
		lines[n + 1] = part2;
	}

	for(n = 0; n < lines.length; n++)
	{
		line = lines[n].trim();
		p = line.lastIndexOf("\x20");
		
		if(p === -1)
		{
			p = line.lastIndexOf("\x09");
			if(p === -1) continue;
		}
		
		part1 = line.substring(0, p + 1).trim();
		part2 = line.substring(p + 1).trim();
		
		num1 = Number(part1);
		num2 = Number(part2);
		
		if(part1.length === 0 || part2.length === 0 || Number.isNaN(num1) || Number.isNaN(num2)) continue;
		
		datapoints[datapoints.length] = {x: num1, y: num2};
		
		// datastr += num1.toExponential() + "\x09" + num2.toExponential() + "\n";
	}
	
	// databox.value = datastr;
	
	delete_all_children(document.getElementById("svg_datapointgroup"));
	
	if(do_update === undefined || do_update)
	{
		// calculate_lsf();
		update(false);
	}
}

function extract_px(str)
{
	return Number(str.substring(0, str.search("px")));
}

function zeropad(number, digits)
{
	var str = number.toFixed(0);
	return "0".repeat(str.length < digits ? digits - str.length : 0) + str;
}

function export_svg()
{
	var figure_div = document.getElementById("figure_div");
	
	var date = new Date(Date.now());
	var filename = "graphic_" + date.getFullYear() + "-" + zeropad(date.getMonth() + 1, 2) + "-" + zeropad(date.getDate(), 2) + "_" + 
			zeropad(date.getHours(), 2) + "-" + zeropad(date.getMinutes(), 2) + "-" + zeropad(date.getSeconds(), 2) + ".svg";
	
	var str = figure_div.innerHTML.replace("width=\"96%\" ", "").replace("height=\"100%\" ", "").replace(/\t/g, "").replace(/\x3C\x2Fpolyline\x3E\x3Cpolyline/g, "</polyline>\n\t<polyline");
	var blob = new File(["<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\" ?>\n", str], filename, {type: "application/svg+xml"});
	var url = URL.createObjectURL(blob);
	window.open(url, "_blank");
}