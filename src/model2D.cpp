#include "model2D.h"

template<> bool Model<char>::flows_to(ArrayCoordinate c1, ArrayCoordinate c2) {
	return ( ( c1.row + directions[this->get(c1.row,c1.col)].row == c2.row ) &&
		 ( c1.col + directions[this->get(c1.row,c1.col)].col == c2.col ) );
}
