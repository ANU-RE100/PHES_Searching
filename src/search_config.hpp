#ifndef SEARCH_CONFIG_H
#define SEARCH_CONFIG_H

#include <string>
#include <iostream>
#include "coordinates.h"


class SearchType {
  public:
    enum type {GREENFIELD, OCEAN, SINGLE_EXISTING, BULK_EXISTING, BULK_PIT};
    constexpr SearchType(type search_type) : value(search_type){}
    constexpr operator type() const { return value; }

    bool existing(){
      return value == SINGLE_EXISTING || value == BULK_EXISTING || value == BULK_PIT;
    }
    bool not_existing(){
      return value == GREENFIELD || value == OCEAN;
    }
    bool grid_cell(){
      return value == GREENFIELD || value == OCEAN || value == BULK_EXISTING || value == BULK_PIT;
    }
    bool single(){
      return value == SINGLE_EXISTING;
    }

    // Filename prefix
    std::string prefix(){
      switch(value){
        case OCEAN:
          return "ocean_";
        case BULK_PIT:
          return "pit_";
        case BULK_EXISTING:
          return "existing_";
        default:
          return "";
      }
    }

    // Used when searching for 8 neighbouring input cells in pairing (i.e. when
    // doing an ocean search we want to read ocean lowers, in all other cases
    // regular neighbours)
    std::string lowers_prefix() {
      switch (value) {
      case OCEAN:
        return "ocean_";
      default:
        return "";
      }
    }

  private:
    type value;
};

class Logger {
  public:
    enum level {DEBUG, ERROR};
    constexpr Logger(level logging_level) : logging_level(logging_level){}
    Logger(char* c){
      if (atoi(c))
        logging_level = DEBUG;
      else
        logging_level = ERROR;
    }
    constexpr operator level() const { return logging_level; }

    void error(std::string message){
      std::cout << message << std::endl;
    }

    bool output_debug(){
      return logging_level == DEBUG;
    }

    void debug(std::string message){
      if (this->output_debug())
        std::cout << message << std::endl;
    }
    void warning(std::string message){
      if (this->output_debug())
        std::cout << message << std::endl;
    }

  private:
    level logging_level = DEBUG;
};

class SearchConfig {
  public:
    SearchType search_type;
    GridSquare grid_square;
    std::string name;
    Logger logger;

    SearchConfig() : search_type(SearchType::GREENFIELD), logger(Logger::ERROR){}
    SearchConfig(int nargs, char **argv) : search_type(SearchType::GREENFIELD), logger(Logger::ERROR) {
      std::string arg1(argv[1]);
      int adj = 0;
      if (arg1.compare("ocean") == 0) {
        search_type = SearchType::OCEAN;
        adj = 1;
        arg1 = argv[1 + adj];
      }
      if (arg1.compare("bulk_existing") == 0) {
        search_type = SearchType::BULK_EXISTING;
        adj = 1;
        arg1 = argv[1 + adj];
      } else if (arg1.compare("bulk_pit") == 0) {
        search_type = SearchType::BULK_PIT;
        adj = 1;
        arg1 = argv[1 + adj];
      } 
      if (arg1.compare("reservoir") == 0) {
        search_type = SearchType::SINGLE_EXISTING;
        adj = 1;
        arg1 = argv[1 + adj];
        if (nargs > 2 + adj)
          logger = Logger(argv[2 + adj]);
      } else {
        try {
          int lon = stoi(arg1);
          grid_square = GridSquare_init(atoi(argv[2 + adj]), lon);
          if (nargs > 3 + adj)
            logger = Logger(argv[3 + adj]);
        } catch (exception &e) {
          search_type = SearchType::SINGLE_EXISTING;
          if (nargs > 2 + adj)
            logger = Logger(argv[2 + adj]);
        }
      }
    }

    std::string filename(){
      if(search_type.grid_cell())
        return search_type.prefix() + str(grid_square);
      return search_type.prefix() + name;
    }
};

extern SearchConfig search_config;

#endif
