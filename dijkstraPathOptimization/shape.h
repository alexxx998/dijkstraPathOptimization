#pragma once
#include "common.h"
#include "edge.h"
#include "vector_math.h"
#include "vertex.h"
#include "face.h"

#include <array>
#include <memory>
#include <unordered_map>
#include <vector>

using namespace std;

class Shape {
public:
	// Set triangulated shape using list of coordinates, list of edges and list of trinagles
	// listOfPoints - list of point coordinates
	// listOfEdges - list of pair of indexes from listOfPoints
	// listOfTriangles - list of triplets of indexes from listOfEdges
	void setShape(const vector<Vector3d>& listOfPoints, const vector<array<int, 2>>& listOfEdges, const vector<array<int, 3>>& listOfTriangles);

	// Dijkstra algorithm
	vector<shared_ptr<Vertex>> findShortestPathByEdges(shared_ptr<Vertex> start, shared_ptr<Vertex> finish) const;

	// Path optimization
	vector<Vector3d> findShortestPath(const vector<shared_ptr<Vertex>>& pathByEdge) const;
	
	shared_ptr<Vertex> getNearestVertex(const Vector3d& point) const;

private:
	vector<shared_ptr<Face>> faces;
};