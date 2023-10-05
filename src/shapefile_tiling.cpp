#include "model2D.h"
#include "phes_base.h"
#include <shapefil.h>

// Reads a list of cells to process from the tasks_file (Eg. 148 -36)
vector<GridSquare> read_tasklist(char *tasks_file) {
  ifstream fd(tasks_file);
  if (!fd) {
    fprintf(stderr, "failed to open task file %s: %s\n", tasks_file, strerror(errno));
    exit(1);
  }
  vector<GridSquare> tasklist;
  string line;
  while (getline(fd, line)) {
    line.erase(remove(line.begin(), line.end(), '\n'), line.end());
    line.erase(remove(line.begin(), line.end(), '\r'), line.end());
    istringstream is(line);
    int lon, lat;
    if (is >> lon)
      if (is >> lat)
        tasklist.push_back(GridSquare_init(lat, lon));
  }
  printf("read %zu tasks\n", tasklist.size());
  fd.close();
  return tasklist;
}

int SHAPE_TYPE = SHPT_POLYGON;

struct Shape {
  vector<GeographicCoordinate> points;
  explicit Shape(vector<GeographicCoordinate>& points) : points(points) {};
  double volume=0;
  string name="";
  int elevation = 0;
};

std::vector<Shape> read_shapefile(string filename, Model<bool> *to_keep, Model<vector<int>> *relevant_polygons, std::string type, int &SHAPE_TYPE){
	vector<Shape> shapes;
  	int iterator = 0;

	SHPHandle SHP = SHPOpen(convert_string(file_storage_location + filename), "rb");
    DBFHandle DBF = DBFOpen(convert_string(file_storage_location + filename), "rb");
    int dbf_field = 0;
    int dbf_name_field = 0;
    int dbf_elevation_field = 0;
    if (type == "RIVER")
      dbf_field = DBFGetFieldIndex(DBF, string("DIS_AV_CMS").c_str());
    if (type == "BLUEFIELD"){
      dbf_field = DBFGetFieldIndex(DBF, string("Vol_total").c_str());
      dbf_elevation_field = DBFGetFieldIndex(DBF, string("Elevation").c_str());
      dbf_name_field = DBFGetFieldIndex(DBF, string("Lake_name").c_str());
    }
    if (SHP != NULL) {
      int nEntities;
      SHPGetInfo(SHP, &nEntities, NULL, NULL, NULL);
      for (int i = 0; i < nEntities; i++) {
        double flow = 0, volume=0;
        if (type == "RIVER"){
          flow = DBFReadDoubleAttribute(DBF, i, dbf_field);
          if (flow < 1)
            continue;
        }
        if (type == "BLUEFIELD"){
          volume = DBFReadDoubleAttribute(DBF, i, dbf_field);
          if (volume < 1)
            continue;
        }
        SHPObject *SHP_shape;
        SHP_shape = SHPReadObject(SHP, i);
        if (SHP_shape == NULL) {
          fprintf(stderr, "Unable to read shape %d, terminating object reading.\n", i);
          break;
        }
        vector<GeographicCoordinate> temp_poly;

        for (int j = 0, iPart = 1; j < SHP_shape->nVertices; j++) {
          if (iPart < SHP_shape->nParts && SHP_shape->panPartStart[iPart] == j) {

            for (int lat = -90; lat < 90; lat++)
              for (int lon = -180; lon < 180; lon++)
                if (to_keep->get(lat + 90, lon + 180)) {
                  relevant_polygons->get_pointer(lat + 90, lon + 180)->push_back(iterator);
                  to_keep->set(lat + 90, lon + 180, false);
                }
            Shape shape = Shape(temp_poly);
            if (type == "RIVER")
              shape.volume = flow;
            if(type=="BLUEFIELD"){
              shape.volume = volume;
              shape.elevation =DBFReadIntegerAttribute(DBF, i, dbf_elevation_field);
              shape.name = string(DBFReadStringAttribute(DBF, i, dbf_name_field));
            }
            shapes.push_back(shape);
            iterator++;
            temp_poly.clear();
            if(type=="BLUEFIELD")
              break;
            iPart++;
          }
          GeographicCoordinate temp_point =
              GeographicCoordinate_init(SHP_shape->padfY[j], SHP_shape->padfX[j]);
          to_keep->set((int)(temp_point.lat + 90), (int)(temp_point.lon + 180), true);
          temp_poly.push_back(temp_point);
        }
        for (int lat = -90; lat < 90; lat++)
          for (int lon = -180; lon < 180; lon++)
            if (to_keep->get(lat + 90, lon + 180)) {
              relevant_polygons->get_pointer(lat + 90, lon + 180)->push_back(iterator);
              to_keep->set(lat + 90, lon + 180, false);
            }
        Shape shape = Shape(temp_poly);
        if (type == "RIVER")
          shape.volume = flow;
        if(type=="BLUEFIELD"){
          shape.volume = volume;
          shape.elevation =DBFReadIntegerAttribute(DBF, i, dbf_elevation_field);
          shape.name = string(DBFReadStringAttribute(DBF, i, dbf_name_field));
        }
        shapes.push_back(shape);
        iterator++;
        SHAPE_TYPE = SHP_shape->nSHPType,
        SHPDestroyObject(SHP_shape);
      }
      printf("%d Polygons imported from %s\n", iterator, convert_string(filename));
    }
    SHPClose(SHP);
	DBFClose(DBF);

	return shapes;
}

void DBFCopy(DBFHandle &oldDBFHandle, DBFHandle &newDBFHandle) {
    int numFields = DBFGetFieldCount(oldDBFHandle);
    int numRecords = DBFGetRecordCount(oldDBFHandle);

    for (int i = 0; i < numFields; i++) {
        char fieldName[12];
        int fieldWidth, fieldDecimals;
        DBFFieldType fieldType = DBFGetFieldInfo(oldDBFHandle, i, fieldName, &fieldWidth, &fieldDecimals);
        DBFAddField(newDBFHandle, fieldName, fieldType, fieldWidth, fieldDecimals);
    }

    for (int i = 0; i < numRecords; i++) {
        for (int j = 0; j < numFields; j++) {
            if (DBFIsAttributeNULL(oldDBFHandle, i, j)) {
                DBFWriteNULLAttribute(newDBFHandle, i, j);
            } else if (DBFGetFieldInfo(oldDBFHandle, j, nullptr, nullptr, nullptr) == FTString) {
                const char* value = DBFReadStringAttribute(oldDBFHandle, i, j);
                DBFWriteStringAttribute(newDBFHandle, i, j, value);
            } else if (DBFGetFieldInfo(oldDBFHandle, j, nullptr, nullptr, nullptr) == FTInteger) {
                int value = DBFReadIntegerAttribute(oldDBFHandle, i, j);
                DBFWriteIntegerAttribute(newDBFHandle, i, j, value);
            } else if (DBFGetFieldInfo(oldDBFHandle, j, nullptr, nullptr, nullptr) == FTDouble) {
                double value = DBFReadDoubleAttribute(oldDBFHandle, i, j);
                DBFWriteDoubleAttribute(newDBFHandle, i, j, value);
            } else if (DBFGetFieldInfo(oldDBFHandle, j, nullptr, nullptr, nullptr) == FTLogical) {
                const char* value = DBFReadLogicalAttribute(oldDBFHandle, i, j);
                DBFWriteLogicalAttribute(newDBFHandle, i, j, *value);
            }
        }
    }
	return;
}


int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cout << "Please specify FILTER, RIVER, or BLUEFIELD." << std::endl;
    exit(1);
  }
  std::string type = std::string(argv[1]);
  if (type != "FILTER" && type != "BLUEFIELD" && type != "RIVER") {
    std::cout << "Please specify FILTER, RIVER, or BLUEFIELD." << std::endl;
    exit(1);
  }
  std::cout << "Tiling started for " << type << std::endl;
  printf("Tiling started\n");

  parse_variables(convert_string("storage_location"));
  parse_variables(convert_string(file_storage_location + "variables"));
  unsigned long start_usec = walltime_usec();

  vector<GridSquare> tasklist = read_tasklist(convert_string(file_storage_location + tasks_file));

  vector<string> shapefile_names;
  if (type == "FILTER")
    shapefile_names = filter_filenames_to_tile;
  else
    shapefile_names.push_back("input/existing_reservoirs/" + existing_reservoirs_shp);    

  string output_location = file_storage_location + "input/shapefile_tiles";
  if(type=="RIVER")
    output_location = file_storage_location + "input/river_shapefile_tiles";
  if(type=="BLUEFIELD")
    output_location = file_storage_location + "input/bluefield_shapefile_tiles";
  mkdir(convert_string(output_location), 0777);
  int z = 0;
  bool first_file = true;

  for (string filename : shapefile_names) {
	
	  Model<vector<int>> *relevant_polygons = new Model<vector<int>>(180, 360);
  	Model<bool> *to_keep = new Model<bool>(180, 360, MODEL_SET_ZERO);
	  std::vector<Shape> shapes = read_shapefile(filename,to_keep,relevant_polygons,type,SHAPE_TYPE);
	
    for (GridSquare gs : tasklist) {
      z++;
      int lat = gs.lat;
      int lon = gs.lon;
      
      string shpName = output_location + "/" + str(gs) + "_shapefile_tile.shp";
      string dbfName = output_location + "/" + str(gs) + "_shapefile_tile.dbf";
      SHPHandle SHP;
      DBFHandle DBF;

      if(first_file){
        SHP = SHPCreate(convert_string(shpName), SHAPE_TYPE);
        DBF = DBFCreate(convert_string(dbfName));
      }
      else{
        SHP = SHPOpen(convert_string(shpName), "rb+");
        DBF = DBFCreate(convert_string(output_location + "/" + str(gs) + "_temp_shapefile_tile.dbf"));
      }		
      
      int field_num=0, elevation_field_num=0, name_field_num =0;
      if(type=="RIVER")
        field_num = DBFAddField(DBF, convert_string("DIS_AV_CMS"), FTDouble, 10, 3);
      if(type=="BLUEFIELD"){
        field_num = DBFAddField(DBF, convert_string("Vol_total"), FTDouble, 10, 3);
        elevation_field_num = DBFAddField(DBF, convert_string("Elevation"), FTInteger, 10, 0);
        name_field_num = DBFAddField(DBF, convert_string("Lake_name"), FTString, 64, 0);
      }
      if (SHP == NULL) {
        printf("Unable to create:%s\n",
            convert_string(output_location + "/" + str(gs) + "_shapefile_tile.shp"));
        exit(1);
      }
      
      for (uint ipolygon = 0; ipolygon < relevant_polygons->get(lat + 90, lon + 180).size();
        ipolygon++) {
        Shape shape = shapes[relevant_polygons->get(lat + 90, lon + 180)[ipolygon]];
        int nVertices = shape.points.size();
        int panParts[1] = {0};
        
        double *padfX = new double[nVertices];
        double *padfY = new double[nVertices];
        double *padfZ = NULL, *padfM = NULL;

        for (int i = 0; i < nVertices; i++) {
          GeographicCoordinate temp = shape.points[i];
          padfX[i] = temp.lon;
          padfY[i] = temp.lat;
        }
        SHPObject *psObject;
        psObject = SHPCreateObject(SHAPE_TYPE, -1, 0, panParts, NULL, nVertices, padfX, padfY,
                      padfZ, padfM);

        int shp_num = SHPWriteObject(SHP, -1, psObject);
        if(type=="RIVER")
          DBFWriteDoubleAttribute(DBF, shp_num, field_num, shape.volume);
        if(type=="BLUEFIELD"){
          DBFWriteDoubleAttribute(DBF, shp_num, field_num, shape.volume);
          DBFWriteIntegerAttribute(DBF, shp_num, elevation_field_num, shape.elevation);
          DBFWriteStringAttribute(DBF, shp_num, name_field_num, shape.name.c_str());
        }
        SHPDestroyObject(psObject);
        delete[] padfY;
        delete[] padfX;
      }
      
      if(!first_file){
        DBFHandle oldDBFHandle = DBFOpen(convert_string(dbfName),"rb");
        DBFCopy(oldDBFHandle, DBF);
        // Delete old DBF file and rename new DBF file to old DBF file
        DBFClose(oldDBFHandle);
        remove(convert_string(dbfName));
        rename(convert_string(output_location + "/" + str(gs) + "_temp_shapefile_tile.dbf"), convert_string(dbfName));
      }

      SHPClose(SHP);
      DBFClose(DBF);
    }
    first_file = false;
    delete to_keep;
    delete relevant_polygons;

  }  

  printf("Tiling finished. Runtime: %.2f sec\n", 1.0e-6 * (walltime_usec() - start_usec));
}