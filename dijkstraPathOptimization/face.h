#pragma once

#include "common.h"
#include "edge.h"
#include "vector_math.h"
#include "vertex.h"

#include <array>
#include <algorithm>
#include <memory>

using namespace std;

struct Face {
	Face(shared_ptr<Edge> e1, shared_ptr<Edge> e2, shared_ptr<Edge> e3) : edges({ e1, e2, e3 }) {}
	Face(array<shared_ptr<Edge>, 3> edges) : edges(edges) {}
	array<shared_ptr<Edge>, 3> edges;

	pair<shared_ptr<Edge>, shared_ptr<Edge>> getAdjacentEdges(shared_ptr<Vertex> vert);
};

