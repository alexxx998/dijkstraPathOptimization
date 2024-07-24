#include "shape.h"

#include <algorithm>
#include <iterator>
#include <unordered_map>
#include <unordered_set>
#include <queue>

using namespace std;

namespace {
    const double PI = 2. * atan2(1, 0);
    const double EPS = 1E-9;
}


void Shape::setShape(const vector<Vector3d>& listOfPoints,
    const vector<array<int, 2>>& listOfEdges,
    const vector<array<int, 3>>& listOfTriangles) {

    vector<shared_ptr<Vertex>> vertices;
    for (const auto& point : listOfPoints) {
        vertices.push_back(make_shared<Vertex>(point));
    }

    vector<shared_ptr<Edge>> edges;
    for (const auto& edge_indices : listOfEdges) {
        auto from = vertices[edge_indices[0]];
        auto to = vertices[edge_indices[1]];
        edges.push_back(make_shared<Edge>(from, to));
    }


    unordered_map<shared_ptr<Vertex>, vector<shared_ptr<Face>>> adjacentFaces;
    for (const auto& triangleId : listOfTriangles) {
        auto edge1 = edges[triangleId[0]];
        auto edge2 = edges[triangleId[1]];
        auto edge3 = edges[triangleId[2]];
        auto face = make_shared<Face>(edge1, edge2, edge3);
        faces.push_back(face);
        unordered_set<shared_ptr<Vertex>> s;
        s.insert(edge1->to);
        s.insert(edge2->to);
        s.insert(edge3->to);
        s.insert(edge1->from);
        s.insert(edge2->from);
        s.insert(edge3->from);
        for (const auto& vert: s)
            adjacentFaces[vert].push_back(face);
    }
    
    for (auto& item : adjacentFaces) {
        auto& curFaces = item.second;
        auto lastEdges = curFaces[0]->getAdjacentEdges(item.first);
        item.first->edges.push_back(lastEdges.first);
        item.first->edges.push_back(lastEdges.second);
        for (int i = 0; i < curFaces.size() - 1; i++) {
            lastEdges = curFaces[i]->getAdjacentEdges(item.first);
            
            unordered_set<shared_ptr<Edge>> used = { lastEdges.first, lastEdges.second};
            
            for (int j = i + 1; j < curFaces.size(); j++) {
                auto curEdges = curFaces[j]->getAdjacentEdges(item.first);
                if (lastEdges.first == curEdges.first || lastEdges.first == curEdges.second
                    || lastEdges.second == curEdges.first || lastEdges.second == curEdges.second){
                    swap(curFaces[i + 1], curFaces[j]);
                    if (used.count(curEdges.first)) {
                        item.first->edges.push_back(curEdges.second);
                        used.insert(curEdges.second);
                    }
                    else {
                        item.first->edges.push_back(curEdges.first);
                        used.insert(curEdges.first);
                    }
                    break;
                }
            }
        }
        item.first->edges.pop_back();
    }
}


vector<shared_ptr<Vertex>> Shape::findShortestPathByEdges(
    shared_ptr<Vertex> start, 
    shared_ptr<Vertex> finish) const{
    unordered_map<shared_ptr<Vertex>, double> distances;
    unordered_map<shared_ptr<Vertex>, shared_ptr<Vertex>> prevVertex;
    priority_queue<pair<double, shared_ptr<Vertex>>,
        vector<pair<double, shared_ptr<Vertex>>>,
        greater<pair<double, shared_ptr<Vertex>>>> queue;

    distances[start] = 0;
    prevVertex[start] = nullptr;
    queue.push({ 0, start });

    while (!queue.empty()) {
        double currentDistance;
        shared_ptr<Vertex> currentVertex;
        tie(currentDistance, currentVertex) = queue.top();
        queue.pop();

        if (currentVertex == finish) {
            break;
        }

        for (const auto& edge : currentVertex->edges) {
            auto neighbor = (currentVertex == edge->from) ? edge->to : edge->from;
            double newDistance = currentDistance + edge->len();

            if (distances.find(neighbor) == distances.end() || newDistance < distances[neighbor]) {
                distances[neighbor] = newDistance;
                prevVertex[neighbor] = currentVertex;
                queue.push({ newDistance, neighbor });
            }
        }
    }

    if (distances.find(finish) == distances.end()) {
        return {};
    }

    vector<shared_ptr<Vertex>> path;
    auto current = finish;
    while (current != nullptr) {
        path.push_back(current);
        current = prevVertex[current];
    }

    reverse(path.begin(), path.end());
    return path;
}

shared_ptr<Vertex> Shape::getNearestVertex(const Vector3d& point) const {
    shared_ptr<Vertex> nearestVertex = nullptr;
    double minDistance = numeric_limits<double>::max();
    for (const auto& face : faces) {
        for (const auto& edge : face->edges) {
            for (const auto& vertex : { edge->from, edge->to }) {
                double distance = (point - vertex->point).len();

                if (distance < minDistance) {
                    minDistance = distance;
                    nearestVertex = vertex;
                }
            }
        }
    }
    return nearestVertex;
}

static bool isSameEdge(shared_ptr<Edge> edge, shared_ptr<Vertex> vert1, shared_ptr<Vertex> vert2) {
    return edge->from == vert1 && edge->to == vert2 || edge->from == vert2 && edge->to == vert1;
}

static Vector3d chooseOutsideDir(shared_ptr<Edge> edge, shared_ptr<Vertex> cornerVert) {
    if (edge->to == cornerVert)
        return (edge->from->point - cornerVert->point);
    return (edge->to->point - cornerVert->point);
}

static double angleBetweenEdges(shared_ptr<Edge> edge1, shared_ptr<Edge> edge2, shared_ptr<Vertex> cornerVert) {
    Vector3d v1 = chooseOutsideDir(edge1, cornerVert), v2 = chooseOutsideDir(edge2, cornerVert);

    return acos((v1 ^ v2) / v1.len() / v2.len());
}

static bool intersectSegment(Vector2d a, Vector2d b, Vector2d c, Vector2d d, Vector2d& left, Vector2d& right) {
    auto intersect_1d = [](double a, double b, double c, double d) {
            if (a > b)  swap(a, b);
            if (c > d)  swap(c, d);
            return max(a, c) <= min(b, d) + EPS;
        };
    auto isBetween = [](double l, double r, double x) {
        return min(l, r) <= x + EPS && x <= max(l, r) + EPS;
        };

    if (!intersect_1d(a.x(), b.x(), c.x(), d.x()) || !intersect_1d(a.y(), b.y(), c.y(), d.y()))
        return false;
    
    struct line {
        double a, b, c;

        line() {}
        line(Vector2d p, Vector2d q) {
            a = p.y() - q.y();
            b = q.x() - p.x();
            c = -a * p.x() - b * p.y();
            norm();
        }

        void norm() {
            double z = sqrt(a * a + b * b);
            if (abs(z) > EPS)
                a /= z, b /= z, c /= z;
        }

        double dist(Vector2d p) const {
            return a * p.x() + b * p.y() + c;
        }
    } m(a, b), n(c, d);
    

    double zn = m.a * n.b - m.b * n.a;
    if (abs(zn) < EPS) {
        if (abs(m.dist(c)) > EPS || abs(n.dist(a)) > EPS)
            return false;
        if (b < a)  swap(a, b);
        if (d < c)  swap(c, d);
        left = max(a, c);
        right = min(b, d);
        return true;
    }
    else {
        left.x() = right.x() = -(m.c * n.b - m.b * n.c) / zn;
        left.y() = right.y() = -(m.a * n.c - m.c * n.a) / zn;
        return isBetween(a.x(), b.x(), left.x())
            && isBetween(a.y(), b.y(), left.y())
            && isBetween(c.x(), d.x(), left.x())
            && isBetween(c.y(), d.y(), left.y());
    }
}


Vector2d getIntersection(const Vector2d& point, const Vector2d& direction, const Vector2d& end) {
    Vector2d left, right;
    auto res = intersectSegment(point, point + direction, Vector2d(0., 0.), end, left, right);
    if (res)
        return left;
    return end;
}

static vector<Vector3d> getShortestSegment(const vector<shared_ptr<Edge>>& edges, shared_ptr<Vertex> cornerVert) {
    vector<Vector2d> projectedVectors;
    projectedVectors.reserve(edges.size());
    projectedVectors.emplace_back(edges[0]->len(), 0.0);
    
    double ang = 0;
    for (int i = 1; i < edges.size(); i++) {
        ang += angleBetweenEdges(edges[i - 1], edges[i], cornerVert);
        projectedVectors.push_back(Vector2d(edges[i]->len(), 0.0).rotate(ang));
    }

    Vector2d point = projectedVectors[0];
    Vector2d dir = projectedVectors.back() - point;
    vector<Vector3d> res;
    double length = 0;
    for (int i = 1; i < edges.size() - 1; i++){
        Vector2d intersection = getIntersection(point, dir, projectedVectors[i]);
        double k = intersection.len() / projectedVectors[i].len();
        res.push_back((chooseOutsideDir(edges[i], cornerVert) * k) + cornerVert->point);
    }
    return res;
}


vector<Vector3d> Shape::findShortestPath(const vector<shared_ptr<Vertex>>& pathByEdge) const {
    vector<Vector3d> res;
    res.push_back(pathByEdge.begin()->get()->point);
    
    for (int i = 1; i+1 < pathByEdge.size(); i += 2) {    
        int firstEdgeId = -1, lastEdgeId = -1;
        for (int j = 0; j < pathByEdge[i]->edges.size(); j++)
        {
            if (isSameEdge(pathByEdge[i]->edges[j], pathByEdge[i], pathByEdge[i - 1]))
                firstEdgeId = j;

            if (isSameEdge(pathByEdge[i]->edges[j], pathByEdge[i], pathByEdge[i + 1]))
                lastEdgeId = j;
        }

        double ang = 0;
        for(int j = min(firstEdgeId, lastEdgeId); j < max(firstEdgeId, lastEdgeId); j++)
            ang += angleBetweenEdges(pathByEdge[i]->edges[j], pathByEdge[i]->edges[j+1], pathByEdge[i]);

        vector<shared_ptr<Edge>> leftEdges, rightEdges;
        for (int j = min(firstEdgeId, lastEdgeId); j <= max(firstEdgeId, lastEdgeId); j++)
            leftEdges.push_back(pathByEdge[i]->edges[j]);
        if (firstEdgeId > lastEdgeId)
            reverse(leftEdges.begin(), leftEdges.end());

        if (firstEdgeId < lastEdgeId) {
            for(int j = firstEdgeId; j >= 0; j--)
                rightEdges.push_back(pathByEdge[i]->edges[j]);
            for(int j = pathByEdge[i]->edges.size() - 1; j >= lastEdgeId; j--)
                rightEdges.push_back(pathByEdge[i]->edges[j]);
        }
        else {
            for (int j = firstEdgeId; j < pathByEdge[i]->edges.size(); j++)
                rightEdges.push_back(pathByEdge[i]->edges[j]);
            for (int j = 0; j <= lastEdgeId; j++)
                rightEdges.push_back(pathByEdge[i]->edges[j]);
        }
        
        auto middlePointsLeft = getShortestSegment(rightEdges, pathByEdge[i]);
        double leftLen = (res.back() - middlePointsLeft[0]).len() + 
            (pathByEdge[i + 1]->point - middlePointsLeft.back()).len();
        for (int j = 1; j < middlePointsLeft.size(); j++)
            leftLen += (middlePointsLeft[j] - middlePointsLeft[j - 1]).len();

        auto middlePointsRight = getShortestSegment(leftEdges, pathByEdge[i]);
        double rightLen = (res.back() - middlePointsRight[0]).len() + 
            (pathByEdge[i + 1]->point - middlePointsRight.back()).len();
        for (int j = 1; j < middlePointsRight.size(); j++)
            rightLen += (middlePointsRight[j] - middlePointsRight[j - 1]).len();

        double currentLen = (pathByEdge[i + 1]->point - pathByEdge[i]->point).len() +
            (pathByEdge[i - 1]->point - pathByEdge[i]->point).len();

        if (currentLen < min(leftLen, rightLen))
            res.push_back(pathByEdge[i]->point);
        else if(leftLen < rightLen)
            copy(middlePointsLeft.begin(), middlePointsLeft.end(), back_inserter(res));
        else
            copy(middlePointsRight.begin(), middlePointsRight.end(), back_inserter(res));

        res.push_back(pathByEdge[i + 1]->point);
    }
    if (pathByEdge.size() % 2 == 0)
        res.push_back(pathByEdge.back()->point);

    return res;
}