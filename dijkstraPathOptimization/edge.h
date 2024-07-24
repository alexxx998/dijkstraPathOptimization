#pragma once

#include "common.h"
#include "face.h"
#include "vector_math.h"
#include "vertex.h"

#include <memory>

using namespace std;

struct Edge {

	Edge(shared_ptr<Vertex> from, shared_ptr<Vertex> to) : from(from), to(to) {};
	double len();
	double lenSqr();

	shared_ptr<Vertex> from, to;
};
	