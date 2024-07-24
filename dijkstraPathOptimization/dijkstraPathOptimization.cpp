#include "shape.h"

#include <iostream>
#include <vector>

using namespace std;

int main()
{
	vector<Vector3d> listOfPoints = {
		{0, 0, 0},
		{1, 0, 0},
		{1, 1, 0},
		{0, 1, 0},
		{0, 0, 1},
		{1, 0, 1},
		{1, 1, 1},
		{0, 1, 1} 
	};

	vector<array<int, 2>> listOfEdges = {
		{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 2},
		{4, 5}, {5, 6}, {6, 7}, {7, 4}, {4, 6},
		{0, 4}, {1, 5}, {2, 6}, {3, 7},
		{0, 7}, {1, 6}, {0, 5}, {2, 7}
	};

	vector<array<int, 3>> listOfTriangles = {
		{0, 1, 4}, {2, 3, 4},
		{5, 6, 9}, {7, 8, 9},
		{10, 14, 8}, {3, 13, 14},
		{1, 12, 15}, {11, 15, 6},
		{0, 11, 16}, {5, 10, 16},
		{13, 2, 17}, {12, 7, 17}
	};

	Shape cube_shape;
	cube_shape.setShape(listOfPoints, listOfEdges, listOfTriangles);

	Vector3d st = { -1, -1, -1 };
	Vector3d fi = { 2, 2, 2 };

	auto start = cube_shape.getNearestVertex(st);
	auto finish = cube_shape.getNearestVertex(fi);

	auto path = cube_shape.findShortestPathByEdges(start, finish);
	double lengthPath = 0;
	for (int i = 1; i < path.size(); i++)
	{
		auto v = path[i]->point - path[i - 1]->point;
		lengthPath += v.len();
	}
	cout.precision(2);
	cout << fixed << "Dijkstra`s path length: " << lengthPath << "\n";
	for (const auto& vert : path) {
		cout << vert->point.x() << " " << vert->point.y() << " " << vert->point.z() << '\n';
	}
	cout << "\n\n\n";

	auto shertestPath = cube_shape.findShortestPath(path);
	lengthPath = 0;
	for (int i = 1; i < shertestPath.size(); i++)
		lengthPath = (shertestPath[i] - shertestPath[i - 1]).len();

	cout << "Optimized Dijkstra`s path length: " << lengthPath << "\n";
	for (const auto& vert : shertestPath) {
		cout << vert.x() << " " << vert.y() << " " << vert.z() << '\n';
	}
}
