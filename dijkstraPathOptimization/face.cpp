#include "face.h"

pair<shared_ptr<Edge>, shared_ptr<Edge>> Face::getAdjacentEdges(shared_ptr<Vertex> vert)
{
	vector<shared_ptr<Edge>> res;
	res.reserve(2);
	for (int i = 0; i < edges.size(); i++)
		if (edges[i]->from == vert || edges[i]->to == vert)
			res.push_back(edges[i]);
	if(res.size() == 2)
		return pair<shared_ptr<Edge>, shared_ptr<Edge>>(res[0], res[1]);
	return pair<shared_ptr<Edge>, shared_ptr<Edge>>();
}
