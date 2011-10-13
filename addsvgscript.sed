#!/bin/sed -f
#
# addsvgscript.sed - Edit the Dunnart SVG in order to add ECMAScript
#                   for hovertext etc.
#
# Alex Stivala
# May 2008
#
# $Id: addsvgscript.sed 4219 2010-10-07 01:43:36Z alexs $
#


/<svg /{
# delete viewBox, width, and height attributes from svg
s/width="[^"]*"//
s/height="[^"]*"//
s/viewBox="[^"]*"//
# insert the onLoad into the main svg element properties
s/<svg \(.*\)>/<svg \1 onload="initLater()">/
}

# add onclick event for tsLabel and tshtLabel objects
s/class="tsLabel" /class="tsLabel" onclick="handleClickEvent(evt)" /
s/class="tshtLabel" /class="tshtLabel" onclick="handleClickEvent(evt)" /

# add the ECMAScript we need straight after the main svg properties 
/<svg .*>/a\
<script type="text/ecmascript" xlink:href="/pro-origami/Window.js"/>\
<script type="text/ecmascript" xlink:href="/pro-origami/mapApp.js"/>\
<script type="text/ecmascript" xlink:href="/pro-origami/helper_functions.js"/>\
<script type="text/ecmascript" xlink:href="/pro-origami/timer.js"/>\
<script type="text/ecmascript" xlink:href="/pro-origami/button.js"/>\
<script type="text/ecmascript" xlink:href="/pro-origami/checkbox_and_radiobutton.js"/>\
<script type="text/ecmascript" xlink:href="/pro-origami/pohovertext.js"/>\
<script type="text/ecmascript" xlink:href="/pro-origami/poselectsse.js"/>\
<script type="text/ecmascript" xlink:href="/pro-origami/SVGPan.js"/>\
<script type="text/ecmascript" xlink:href="/pro-origami/poqptabsearch.js"/>

# add the SVG group 'cartoon' open immediately after that
/<svg .*>/a\
<g id="cartoon">

# add the cartoon group close and new groups used by mapApp.js at the end
# just before the close svg
/<\/svg>/i\
</g>\
<defs>\
		<!-- symbols for check boxes -->\
		<symbol id="checkBoxRect" overflow="visible">\
			<rect x="-5" y="-5" width="10" height="10" fill="white" stroke="dimgray"\
				stroke-width="1" cursor="pointer"/>\
		</symbol>\
		<symbol id="checkBoxCross" overflow="visible">\
			<g pointer-events="none" stroke="dimgray" stroke-width="1">\
				<line x1="-3" y1="-3" x2="3" y2="3"/>\
				<line x1="3" y1="-3" x2="-3" y2="3"/>\
			</g>\
		</symbol>\
		<!-- symbols for radio buttons -->\
		<symbol id="radioBorder" overflow="visible">\
			<circle fill="white" stroke="dimgray" stroke-width="1.5" r="5"/>\
		</symbol>\
		<symbol id="radioPoint" overflow="visible">\
			<circle fill="dimgray" r="3" pointer-events="none"/>\
		</symbol>\
	</defs>\
<g id="checkboxes" />\
<g id="radioButtonsSelectmode" />\
<g id="textbuttons" />\
<g id="toolTip" />

