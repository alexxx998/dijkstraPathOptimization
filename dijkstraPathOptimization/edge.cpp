#include "edge.h"

double Edge::len() {
	return (this->from->point - this->to->point).len();
}

double Edge::lenSqr() {
	return (this->from->point - this->to->point).lenSqr();
}
