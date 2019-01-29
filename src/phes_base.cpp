#include <sys/time.h>

//Being risky
#include <bits/stdc++.h>
using namespace std;

#include "phes_base.h"
#include "model2D.h"
#include "TIFF_IO.h"

int convert_to_int(double f)
{
	if(f>=0)
		return (int) (f+0.5);
	else
		return (int) (f-0.5);
}

double max_over_wall_heights(double *a)
{
	double amax = -1.0e20;
	for (int ih=0; ih<NWALL_HEIGHTS; ih++)
		amax = MAX(amax, a[ih]);

	return amax;
}

double convert_to_dam_volume(int height, double length)
{
	return (((height+freeboard)*(cwidth+dambatter*(height+freeboard)))/1000000)*length;
}

double convert_to_dam_volume(double height, double length)
{
	return (((height+freeboard)*(cwidth+dambatter*(height+freeboard)))/1000000)*length;
}

double linear_interpolate(double value, const double *x_values, const double *y_values)
{
	int i = 0;
	while (x_values[i]<value) {
		if (i==NWALL_HEIGHTS-1)
			return INF;
		else
			i++;
	}

	double xlower = (i) ? x_values[i-1] : 0;
	double ylower = (i) ? y_values[i-1] : 0;
	double r = x_values[i]-xlower;

	return (ylower+(y_values[i]-ylower)*(value-xlower)/r);
}

string str(int i)
{
	char buf[32];
	sprintf(buf, "%d", i);
	string to_return(buf);
	return to_return;
}

unsigned long walltime_usec()
{
	struct timeval now;
	gettimeofday(&now,(struct timezone*)0);
	return (1000000*now.tv_sec + now.tv_usec);
}

double find_required_volume(int energy, int head)
{
	return (((double)(energy)*J_GWh_conversion)/((double)(head)*water_density*gravity*generation_efficiency*usable_volume*cubic_metres_GL_conversion));
}

char* convert_string(string str){
	char *c = new char[str.length() + 1];
	strcpy(c, str.c_str());
	return c;
}

string ReplaceAll(string str, const string& from, const string& to) {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    return str;
}

void write_to_csv_file(FILE *csv_file, vector<string> cols){
	for(uint i = 0; i<cols.size();i++){
		if (cols[i].find(',') != std::string::npos){
			cols[i] = ReplaceAll(cols[i], string("\""), string("\"\""));
			cols[i] = ReplaceAll(cols[i], string("  "), string(""));
			cols[i] = ReplaceAll(cols[i], string("\n"), string(""));
			cols[i] = '"'+cols[i]+'"';
		}
		fprintf(csv_file, "%s", convert_string(cols[i]));
		if(i!=cols.size()-1)
			fprintf(csv_file, ",");
	}
	fprintf(csv_file, "\n");
}

vector<string> read_from_csv_file(string line){
	vector<string> cols;
	istringstream ss(line);
	string col;
    while (getline(ss, col, ',')) {
        cols.push_back(col);
    }
	return cols;
}

string dtos(double f, int nd) {
	stringstream ss;
	ss << fixed << std::setprecision(nd) << f;
	return ss.str();
}

	
Model_int16* read_DEM_with_borders(GridSquare sc){
	int b[2] = {0,0};
	Model_int16 *DEM = Model_int16_create(b, MODEL_SET_ZERO);
	const int neighbors[9][4][2] = {
		//[(Tile coordinates) , (Tile base)		 		  , (Tile limit)				  , (Tile offset)	 	       ]
		{ {sc.lat  ,sc.lon  } , {border,      border	 }, {border+3600,  	3600+border	 }, {border-1,    border     } },
		{ {sc.lat+1,sc.lon-1} , {0,			  0		 	 }, {border, 	    border	 	 }, {border-3601, border-3600} },
		{ {sc.lat+1,sc.lon  } , {0,	      	  border	 }, {border,	    3600+border	 }, {border-3601, border     } },
		{ {sc.lat+1,sc.lon+1} , {0,	      	  3600+border}, {border,        3600+2*border}, {border-3601, border+3600} },
		{ {sc.lat  ,sc.lon+1} , {border-1,    3600+border}, {3600+border,   3600+2*border}, {border-1,    border+3600} },
		{ {sc.lat-1,sc.lon+1} , {3600+border, 3600+border}, {3600+2*border, 3600+2*border}, {border+3599, border+3600} },
		{ {sc.lat-1,sc.lon  } , {3600+border, border	 }, {3600+2*border, 3601+border	 }, {border+3599, border     } },
		{ {sc.lat-1,sc.lon-1} , {3600+border, 0		 	 }, {3600+2*border, border	 	 }, {border+3599, border-3600} },
		{ {sc.lat  ,sc.lon-1} , {border-1,    0		 	 }, {3600+border,   border	 	 }, {border-1,    border-3600} }
	};

	for (int i=0; i<9; i++) {
		GridSquare gs = GridSquare_init(neighbors[i][0][0], neighbors[i][0][1]);
		ArrayCoordinate tile_start = ArrayCoordinate_init(neighbors[i][1][0], neighbors[i][1][1], get_origin(gs, border));
		ArrayCoordinate tile_end = ArrayCoordinate_init(neighbors[i][2][0], neighbors[i][2][1], get_origin(gs, border));
		ArrayCoordinate tile_offset = ArrayCoordinate_init(neighbors[i][3][0], neighbors[i][3][1], get_origin(gs, border));

		try{
			Model_int16 *DEM_temp = TIFF_Read_int16(convert_string("input/"+str(gs)+"_1arc_v3.tif"), NULL, NULL);
			if (i==0) {
				int bordered_shape[2] = {DEM_temp->shape[0]+2*border-1,DEM_temp->shape[1]+2*border-1};
				DEM = Model_int16_create(bordered_shape, MODEL_SET_ZERO);
			}
			for(int row = tile_start.row ; row < tile_end.row ; row++)
				for(int col = tile_start.col ; col < tile_end.col; col++)
					DEM->d[row][col] = DEM_temp->d[row-tile_offset.row][col-tile_offset.col];
			Model_int16_free(DEM_temp);
		}catch (int e){
			if(display)
				fprintf(stderr, "Could not find file %s: %s\n", convert_string("input/"+str(gs)+"_1arc_v3.tif"), strerror(errno));
			if (i==0)
				throw(1);
		}
	}
	return DEM;
}

Models Models_init(GridSquare sc){
	Models models;
	GridSquare neighbors[9] = {
		(GridSquare){sc.lat  ,sc.lon  },
		(GridSquare){sc.lat+1,sc.lon-1},
		(GridSquare){sc.lat+1,sc.lon  },
		(GridSquare){sc.lat+1,sc.lon+1},
		(GridSquare){sc.lat  ,sc.lon+1},
		(GridSquare){sc.lat-1,sc.lon+1},
		(GridSquare){sc.lat-1,sc.lon  },
		(GridSquare){sc.lat-1,sc.lon-1},
		(GridSquare){sc.lat  ,sc.lon-1}};
	for(int i = 0; i<9; i++){
		models.neighbors[i] = neighbors[i];
	}
	models.origin = get_origin(neighbors[1], border);
	return models;
}

void set_FOM(Pair* pair){
	pair->FOM = pair->head * pair->upper.water_rock * pair->lower.water_rock / (pair->upper.water_rock + pair->lower.water_rock) - 125.0/SQRT(pair->slope);
}

string str(Test test){
	return to_string(test.energy_capacity)+"GWh_"+to_string(test.storage_time)+"h";
}

bool file_exists (char* name) {
	ifstream infile(name);
    return infile.good();
}