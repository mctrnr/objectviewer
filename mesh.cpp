#include <iostream>
#include <vector>

#include "./mesh.h"

using namespace std;

Mesh::Mesh() {
  _cur_mtl = -1;
}

// This will be called by the obj parser
void Mesh::AddVertex(const Vec3f& v) {
  Vertex vert;
  vert.v = v;
  vert.n = Vec3f::makeVec(0, 1, 0);
  vertices.push_back(vert);
  // updates the bounding box
  _bb(v);
}

// This will be called by the obj parser
void Mesh::AddTextureVertex(const Vec3f& v) {
  texVertices.push_back(v);
}

// p is the list of indices of vertices for this polygon.  For example,
// if p = {0, 1, 2} then the polygon is a triangle with the zeroth, first and
// second vertices added using AddVertex.
//
// pt is the list of texture indices for this polygon, similar to the
// actual vertices described above.
void Mesh::AddPolygon(const std::vector<int>& p, const std::vector<int>& pt) {
  int vertSize = p.size();
  int texSize = pt.size();
  vector<int> vPolygon;
  vector<int> tPolygon;
  if (vertSize != texSize) {
    cout << "Error: Vertices size doesn't match Texture\n";
  }
  for (int i = 0; i < vertSize; i++) {
    if (vertices.size() < p[i]) cout << "Vert index error\n";
    vertices[p[i]].faces.push_back(polygons.size());
    vPolygon.push_back(p[i]);
    if (texVertices.size() < pt[i]) cout << "texVert index error\n";
    tPolygon.push_back(pt[i]);
  }
  polygons.push_back(vPolygon);
  textPolygons.push_back(tPolygon);
  Vec3f normal = Vec3f::makeVec(0, 1, 0);
  if (vPolygon.size() >=3) {
    Vec3f a = vertices[vPolygon[0]].v - vertices[vPolygon[1]].v;
    Vec3f b = vertices[vPolygon[0]].v - vertices[vPolygon[2]].v;
    normal = a^b;
    normal = normal.unit();
  }
  if (vPolygon.size() >= 5) cout << vPolygon.size() << endl;
  surfaceNormals.push_back(normal);
  // updates the poly2mat map
  _polygon2material.push_back(_cur_mtl);
}

// Computes a normal for each vertex.
void Mesh::compute_normals() {
  Vec3f zero = Vec3f::makeVec(0, 0, 0);
  Vec3f up = Vec3f::makeVec(0, 1, 0);
  for (int i = 0; i < vertices.size(); i++) {
    Vec3f sum = Vec3f::makeVec(0, 0, 0);
    for (int j = 0; j < vertices[i].faces.size(); j++) {
      sum += surfaceNormals[vertices[i].faces[j]];
    }
    sum /= vertices[i].faces.size();
    if (sum == zero) {
      sum = surfaceNormals[vertices[i].faces[0]];
    }
    vertices[i].n = sum.unit();
  }
}
