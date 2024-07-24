#pragma once
#include "common.h"
#include "vector_math.h"
#include "edge.h"

#include <memory>
#include <vector>

using namespace std;

struct Vertex {
	Vertex(Vector3d point) : point(point) {}
	Vertex(double x, double y, double z) : point(x, y, z) {}

	Vector3d point;
	vector<shared_ptr<Edge>> edges;
};