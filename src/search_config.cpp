#include "search_config.hpp"

SearchConfig search_config;

string format_for_filename(string s){
	replace(s.begin(), s.end(), ' ' , '_');
	s.erase(remove(s.begin(), s.end(), '"'), s.end());
	return s;
}
