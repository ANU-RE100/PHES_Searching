#include "kml.h"

string kml_start = 
"<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n"
"  <Document id=\"Layers\">\n"
"    <Style id=\"reservoir_style\">\n"
"      <LineStyle>\n"
"        <width>0</width>\n"
"      </LineStyle>\n"
"    </Style>\n"
"    <Style id=\"wall_style\">\n"
"      <LineStyle>\n"
"        <width>0</width>\n"
"      </LineStyle>\n"
"      <PolyStyle>\n"
"        <color>ff999999</color>\n"
"      </PolyStyle>\n"
"    </Style>\n"
"    <Style id=\"pipe_style\">\n"
"      <LineStyle>\n"
"        <color>88ffffff</color>\n"
"        <width>2</width>\n"
"      </LineStyle>\n"
"    </Style>\n"
"    <Style id=\"pin_style\">\n"
"      <LabelStyle>\n"
"        <color>00ffffff</color>\n"
"        <scale>0.6</scale>\n"
"      </LabelStyle>\n"
"    </Style>\n";

string kml_end = 
"  </Document>\n"
"</kml>\n";

string get_html(Reservoir* reservoir, Pair* pair){
	string newline = "";
	if(reservoir->pit){
		return
"              <html>"+newline+
"              <head><META http-equiv=\"Content-Type\" content=\"text/html\"><meta http-equiv=\"content-type\" content=\"text/html; charset=UTF-8\"></head>"+newline+
"              <body style=\"margin:0px 0px 0px 0px;overflow:auto;background:#FFFFFF;\">"+newline+
"              <table style=\"font-family:Arial,Verdana,Times;font-size:12px;text-align:left;width:100%;border-collapse:collapse;padding:3px 3px 3px 3px\">"+newline+
"              <tr style=\"text-align:center;font-weight:bold;background:#9CBCE2\"><td>"+reservoir->identifier+"</td></tr>"+newline+
"              <tr><td>"+newline+
"              <table style=\"font-family:Arial,Verdana,Times;font-size:12px;text-align:left;width:100%;border-spacing:0px; padding:3px 3px 3px 3px\">"+newline+
(pair != NULL ? "              <tr><td>Class</td><td>"+string(1,pair->category)+"</td></tr>"+newline : "") +
"              <tr bgcolor=\"#D4E4F3\"><td>Elevation (m)</td><td>"+to_string(reservoir->elevation)+"</td></tr>"+newline+
"              <tr><td>Depth (m)</td><td>"+dtos(reservoir->dam_height,0)+"</td></tr>"+newline+
"              <tr bgcolor=\"#D4E4F3\"><td>Latitude</td><td>"+dtos(reservoir->latitude,4)+"</td></tr>"+newline+
"              <tr><td>Longitude</td><td>"+dtos(reservoir->longitude,4)+"</td></tr>"+newline+
"              <tr bgcolor=\"#D4E4F3\"><td>Volume (GL)</td><td>"+dtos(reservoir->volume,1)+"</td></tr>"+newline+
"              <tr><td>Country</td><td>"+reservoir->country+"</td></tr>"+newline+
"              </table></td></tr></table></body></html>\n";
	}else if(reservoir->brownfield){
		return
"              <html>"+newline+
"              <head><META http-equiv=\"Content-Type\" content=\"text/html\"><meta http-equiv=\"content-type\" content=\"text/html; charset=UTF-8\"></head>"+newline+
"              <body style=\"margin:0px 0px 0px 0px;overflow:auto;background:#FFFFFF;\">"+newline+
"              <table style=\"font-family:Arial,Verdana,Times;font-size:12px;text-align:left;width:100%;border-collapse:collapse;padding:3px 3px 3px 3px\">"+newline+
"              <tr style=\"text-align:center;font-weight:bold;background:#9CBCE2\"><td>"+reservoir->identifier+"</td></tr>"+newline+
"              <tr><td>"+newline+
"              <table style=\"font-family:Arial,Verdana,Times;font-size:12px;text-align:left;width:100%;border-spacing:0px; padding:3px 3px 3px 3px\">"+newline+
//"              <tr><td>Reservoir Ref.</td><td>"+reservoir->identifier+"</td></tr>"+newline+
(pair != NULL ? "              <tr><td>Class</td><td>"+string(1,pair->category)+"</td></tr>"+newline : "")+
"              <tr bgcolor=\"#D4E4F3\"><td>Elevation (m)</td><td>"+to_string(reservoir->elevation)+"</td></tr>"+newline+
"              <tr><td>Latitude</td><td>"+dtos(reservoir->latitude,4)+"</td></tr>"+newline+
"              <tr bgcolor=\"#D4E4F3\"><td>Longitude</td><td>"+dtos(reservoir->longitude,4)+"</td></tr>"+newline+
"              <tr><td>Volume (GL)</td><td>"+dtos(reservoir->volume,1)+"</td></tr>"+newline+
"              <tr bgcolor=\"#D4E4F3\"><td>Country</td><td>"+reservoir->country+"</td></tr>"+newline+
"              </table></td></tr></table></body></html>\n";
	}else{
		return
"              <html>"+newline+
"              <head><META http-equiv=\"Content-Type\" content=\"text/html\"><meta http-equiv=\"content-type\" content=\"text/html; charset=UTF-8\"></head>"+newline+
"              <body style=\"margin:0px 0px 0px 0px;overflow:auto;background:#FFFFFF;\">"+newline+
"              <table style=\"font-family:Arial,Verdana,Times;font-size:12px;text-align:left;width:100%;border-collapse:collapse;padding:3px 3px 3px 3px\">"+newline+
"              <tr style=\"text-align:center;font-weight:bold;background:#9CBCE2\"><td>"+reservoir->identifier+"</td></tr>"+newline+
"              <tr><td>"+newline+
"              <table style=\"font-family:Arial,Verdana,Times;font-size:12px;text-align:left;width:100%;border-spacing:0px; padding:3px 3px 3px 3px\">"+newline+
//"              <tr><td>Reservoir Ref.</td><td>"+reservoir->identifier+"</td></tr>"+newline+
(pair != NULL ? "              <tr><td>Class</td><td>"+string(1,pair->category)+"</td></tr>"+newline : "")+
"              <tr bgcolor=\"#D4E4F3\"><td>Elevation (m)</td><td>"+to_string(reservoir->elevation)+"</td></tr>"+newline+
"              <tr><td>Latitude</td><td>"+dtos(reservoir->latitude,4)+"</td></tr>"+newline+
"              <tr bgcolor=\"#D4E4F3\"><td>Longitude</td><td>"+dtos(reservoir->longitude,4)+"</td></tr>"+newline+
"              <tr><td>Area (ha)</td><td>"+dtos(reservoir->area,0)+"</td></tr>"+newline+
"              <tr bgcolor=\"#D4E4F3\"><td>Volume (GL)</td><td>"+dtos(reservoir->volume,1)+"</td></tr>"+newline+
"              <tr><td>Dam Wall Height (m)</td><td>"+dtos(reservoir->dam_height,1)+"</td></tr>"+newline+
"              <tr bgcolor=\"#D4E4F3\"><td>Dam Length (m)</td><td>"+dtos(reservoir->dam_length,0)+"</td></tr>"+newline+
"              <tr><td>Dam Volume (GL)</td><td>"+dtos(reservoir->dam_volume,1)+"</td></tr>"+newline+
"              <tr bgcolor=\"#D4E4F3\"><td>Water/Rock Ratio</td><td>"+dtos(reservoir->water_rock,1)+"</td></tr>"+newline+
"              <tr><td>Country</td><td>"+reservoir->country+"</td></tr>"+newline+
"              </table></td></tr></table></body></html>\n";
	}
}

string get_html(Pair* pair){
	string newline = "";
	return
"              <html>"+newline+
"              <head><META http-equiv=\"Content-Type\" content=\"text/html\"><meta http-equiv=\"content-type\" content=\"text/html; charset=UTF-8\"></head>"+newline+
"              <body style=\"margin:0px 0px 0px 0px;overflow:auto;background:#FFFFFF;\">"+newline+
"              <table style=\"font-family:Arial,Verdana,Times;font-size:12px;text-align:left;width:100%;border-collapse:collapse;padding:3px 3px 3px 3px\">"+newline+
"              <tr style=\"text-align:center;font-weight:bold;background:#9CBCE2\">"+newline+
"              <td>"+pair->identifier+"</td></tr>"+newline+
"              <tr><td>"+newline+
"              <table style=\"font-family:Arial,Verdana,Times;font-size:12px;text-align:left;width:100%;border-spacing:0px; padding:3px 3px 3px 3px\">"+newline+
"              <tr><td>Class</td><td>"+string(1,pair->category)+"</td></tr>"+newline+
"              <tr bgcolor=\"#D4E4F3\"><td>Head (m)</td><td>"+to_string(pair->head)+"</td></tr>"+newline+
"              <tr><td>Separation (km)</td><td>"+dtos(pair->distance,1)+"</td></tr>"+newline+
"              <tr bgcolor=\"#D4E4F3\"><td>Average Slope (%)</td><td>"+dtos(pair->slope*100,0)+"</td></tr>"+newline+
"              <tr><td>Volume (GL)</td><td>"+dtos(pair->volume,1)+"</td></tr>"+newline+
"              <tr bgcolor=\"#D4E4F3\"><td>Water to Rock (Pair)</td><td>"+dtos(pair->water_rock,1)+"</td></tr>"+newline+
"              <tr><td>Energy (GWh)</td><td>"+energy_capacity_to_string(pair->energy_capacity)+"</td></tr>"+newline+
"              <tr bgcolor=\"#D4E4F3\"><td>Storage time (h)</td><td>"+to_string(pair->storage_time)+"</td></tr>"+newline+
"              <tr><td>Country</td><td>"+pair->country+"</td></tr>"+newline+
"              </table>"+newline+
"              </body>"+newline+
"              </html>\n";
}

string get_reservoir_geometry(Reservoir_KML_Coordinates coordinates){
	return
"        <MultiGeometry>\n"
"          <Polygon>\n"
"            <extrude>1</extrude>\n"
"            <altitudeMode>absolute</altitudeMode>\n"
"            <outerBoundaryIs><LinearRing>\n"
"              <coordinates>"+coordinates.reservoir+"</coordinates>\n"
"            </LinearRing></outerBoundaryIs>\n"
"          </Polygon>\n"
"        </MultiGeometry>\n";
}

string get_dam_geometry(Reservoir_KML_Coordinates coordinates){
	string to_return =
"        <MultiGeometry>\n";
	if(coordinates.is_turkeys_nest){
		to_return+=
"          <Polygon>\n"
"            <extrude>1</extrude>\n"
"            <altitudeMode>absolute</altitudeMode>\n"
"            <outerBoundaryIs><LinearRing>\n"
"              <coordinates>"+coordinates.dam[0]+"</coordinates>\n"
"            </LinearRing></outerBoundaryIs>\n"
"            <innerBoundaryIs><LinearRing>\n"
"              <coordinates>"+coordinates.reservoir+"</coordinates>\n"
"            </LinearRing></innerBoundaryIs>\n"
"          </Polygon>\n";
	}else{
		for(uint i = 0; i<coordinates.dam.size(); i++){
			to_return+=
"          <Polygon>\n"
"            <extrude>1</extrude>\n"
"            <altitudeMode>absolute</altitudeMode>\n"
"            <outerBoundaryIs><LinearRing>\n"
"              <coordinates>"+coordinates.dam[i]+"</coordinates>\n"
"            </LinearRing></outerBoundaryIs>\n"
"          </Polygon>\n";
		}
	}
	return to_return+
"        </MultiGeometry>\n";
}

string get_line_geometry(string coordinates){
	return
"        <LineString>\n"
"          <tessellate>1</tessellate>\n"
"          <extrude>1</extrude>\n"
"          <altitudeMode>clampToGround</altitudeMode>\n"
"          <coordinates>"+coordinates+"</coordinates>\n"
"        </LineString>\n";
}

string get_point_geometry(string coordinates){
	return
"        <Point>\n"
"          <coordinates>"+coordinates+"</coordinates>\n"
"        </Point>\n";
}

string get_reservoir_kml(Reservoir* reservoir, string colour, Reservoir_KML_Coordinates coordinates, Pair* pair){
	return
"      <Placemark>\n"
"        <name><![CDATA["+reservoir->identifier+"]]></name>\n"
"        <styleUrl>#reservoir_style</styleUrl>\n"
"        <Style>\n"
"          <PolyStyle>\n"
"            <color>"+colour+"</color>\n"
"          </PolyStyle>\n"
"        </Style>\n"
			+get_reservoir_geometry(coordinates)+
"        <description>\n"
"           <![CDATA[\n"
			+get_html(reservoir, pair)+
"           ]]>\n"
"        </description>\n"
"      </Placemark>\n";
}

string get_dam_kml(Reservoir* reservoir, Reservoir_KML_Coordinates coordinates){
	return
"      <Placemark>\n"
"        <name><![CDATA["+reservoir->identifier+" Dam]]></name>\n"
"        <styleUrl>#wall_style</styleUrl>\n"
			+get_dam_geometry(coordinates)+
"      </Placemark>\n";
}

string get_line_kml(Pair* pair, string coordinates){
	return
"      <Placemark>\n"
"        <name><![CDATA["+pair->identifier+"]]></name>\n"
"        <styleUrl>#pipe_style</styleUrl>\n"
			+get_line_geometry(coordinates)+
"        <description>\n"
"           <![CDATA[\n"
				+get_html(pair)+
"           ]]>\n"
"        </description>\n"
"      </Placemark>\n";
}

string hex(int c){
	stringstream to_return;
	to_return << setw(2) << setfill('0') << hex << c;
	return to_return.str();
}

string get_colour(char category){
	int opacity = convert_to_int(good_colour[0]+(((category-'A')/4.0)*(bad_colour[0]-good_colour[0])));
	int blue = convert_to_int(good_colour[1]+(((category-'A')/4.0)*(bad_colour[1]-good_colour[1])));
	int green = convert_to_int(good_colour[2]+(((category-'A')/4.0)*(bad_colour[2]-good_colour[2])));
	int red = convert_to_int(good_colour[3]+(((category-'A')/4.0)*(bad_colour[3]-good_colour[3])));
	string to_return = hex(opacity)+hex(blue)+hex(green)+hex(red);
	return to_return;
}

string get_scale(char category){
	for(CategoryCutoff category_cutoff:category_cutoffs)
		if(category==category_cutoff.category)
			return dtos(1.25-((category-'A')/8.0),3);
	return "0";
}

string get_point_kml(Pair* pair, string coordinates){
	return
"      <Placemark>\n"
"        <name><![CDATA["+pair->identifier+"]]></name>\n"
"        <styleUrl>#pin_style</styleUrl>\n"
"        <Style>\n"
"          <IconStyle>\n"
"            <scale>"+get_scale(pair->category)+"</scale>\n"
"            <color>"+get_colour(pair->category)+"</color>\n"
"            <Icon><href>http://maps.google.com/mapfiles/kml/paddle/"+pair->category+".png</href></Icon>\n"
"          </IconStyle>\n"
"        </Style>\n"
			+get_point_geometry(coordinates)+
"        <altitudeMode>relativeToGround</altitudeMode>\n"
"        <description>\n"
"           <![CDATA[\n"
				+get_html(pair)+
"           ]]>\n"
"        </description>\n"
"      </Placemark>\n";
}

string generate_folder(string name, vector<string> records){
	string to_return = 
"    <Folder>\n"
"      <name>"+name+"</name>\n";
	for(uint i = 0; i<records.size();i++)
		to_return+=records[i];
	to_return+="    </Folder>\n";
	return to_return;
}

string output_kml(KML_Holder* kml_holder, string square, Test test){
	string to_return;
	to_return+=kml_start;
	to_return+=generate_folder(square+" Upper Reservoirs "+str(test), kml_holder->uppers);
	to_return+=generate_folder(square+" Upper Dams "+str(test), kml_holder->upper_dams);
	to_return+=generate_folder(square+" Lower Reservoirs "+str(test), kml_holder->lowers);
	to_return+=generate_folder(square+" Lower Dams "+str(test), kml_holder->lower_dams);
	to_return+=generate_folder(square+" Lines "+str(test), kml_holder->lines);
	to_return+=generate_folder(square+" Points "+str(test), kml_holder->points);
	to_return+=kml_end;
	return to_return;
}

void update_kml_holder(KML_Holder* kml_holder, Pair* pair, Pair_KML* pair_kml, bool keep_upper, bool keep_lower){
	kml_holder->points.push_back(get_point_kml(pair, pair_kml->point));
	kml_holder->lines.push_back(get_line_kml(pair, pair_kml->line));
	if (keep_upper) {
		kml_holder->uppers.push_back(get_reservoir_kml(&pair->upper, upper_colour, pair_kml->upper, pair));
		if(!pair->upper.brownfield)
			kml_holder->upper_dams.push_back(get_dam_kml(&pair->upper, pair_kml->upper));
	}
	if (keep_lower) {
		if(!pair->lower.ocean)
			kml_holder->lowers.push_back(get_reservoir_kml(&pair->lower, lower_colour, pair_kml->lower, pair));
		if(!pair->lower.brownfield && !pair->lower.ocean)
			kml_holder->lower_dams.push_back(get_dam_kml(&pair->lower, pair_kml->lower));
	}
}
