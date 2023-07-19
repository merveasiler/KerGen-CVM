// @author Merve Asiler

#include "Painter.h"
#include "BaseGeoOpUtils.h"

void Painter::getColorfulShapeSep(Mesh* mesh, SoSeparator* res, vector<double> colorSource)
{
	// Gouraud shading
	SoShapeHints* hints = new SoShapeHints;
	hints->creaseAngle = 3.14;
	res->addChild(hints);

	SoMaterial* mat = new SoMaterial();
	mat->transparency = 0;
	res->addChild(mat);

	SoCoordinate3* coords = new SoCoordinate3();
	for (int c = 0; c < mesh->getNumOfVerts(); c++) {
		coords->point.set1Value(c, mesh->getVertex(c).coords[0], mesh->getVertex(c).coords[1], mesh->getVertex(c).coords[2]);
		mat->diffuseColor.set1Value(c, mesh->getVertex(c).color[0], mesh->getVertex(c).color[1], mesh->getVertex(c).color[2]);
		//mat->diffuseColor.set1Value(c, colorSource[c], 0, 1.0 - colorSource[c]);
	}

	SoMaterialBinding* materialBinding = new SoMaterialBinding; //for 2+ diffuse color usage on the same mesh
	materialBinding->value = SoMaterialBinding::PER_VERTEX_INDEXED;
	res->addChild(materialBinding);

	SoIndexedFaceSet* faceSet = new SoIndexedFaceSet();
	int start_index = mesh->getNumOfVerts();
	for (int c = 0; c < mesh->getNumOfTris(); c++)
	{
		faceSet->coordIndex.set1Value(start_index + c * 4, mesh->getTriangle(c).corners[0]);
		faceSet->coordIndex.set1Value(start_index + c * 4 + 1, mesh->getTriangle(c).corners[1]);
		faceSet->coordIndex.set1Value(start_index + c * 4 + 2, mesh->getTriangle(c).corners[2]);
		faceSet->coordIndex.set1Value(start_index + c * 4 + 3, -1);
	}

	res->addChild(coords);
	res->addChild(faceSet);
}

void Painter::getShapeSep(Mesh* mesh, SoSeparator* res)
{
	// Paint all vertices with the same color
	SoMaterial* mat = new SoMaterial();
	mat->diffuseColor.setValue(1, 1, 1);
	mat->transparency = 0.5;
	res->addChild(mat);

	// Gouraud shading
	SoShapeHints* hints = new SoShapeHints;
	hints->creaseAngle = 3.14;
	res->addChild(hints);

	SoCoordinate3* coords = new SoCoordinate3();
	for (int c = 0; c < mesh->getNumOfVerts(); c++)
		coords->point.set1Value(c, mesh->getVertex(c).coords[0],
			mesh->getVertex(c).coords[1],
			mesh->getVertex(c).coords[2]);

	SoIndexedFaceSet* faceSet = new SoIndexedFaceSet();
	for (int c = 0; c < mesh->getNumOfTris(); c++)
	{
		faceSet->coordIndex.set1Value(c * 4, mesh->getTriangle(c).corners[0]);
		faceSet->coordIndex.set1Value(c * 4 + 1, mesh->getTriangle(c).corners[1]);
		faceSet->coordIndex.set1Value(c * 4 + 2, mesh->getTriangle(c).corners[2]);
		faceSet->coordIndex.set1Value(c * 4 + 3, -1);
	}

	res->addChild(coords);
	res->addChild(faceSet);
}

void Painter::getShapeSep(Mesh* mesh, MaterialSetting* materialSetting, SoSeparator* res)
{
	// Paint all vertices with the same color
	SoMaterial* mat = new SoMaterial();
	mat->diffuseColor.setValue(materialSetting->getRed(), materialSetting->getGreen(), materialSetting->getBlue());		// 0.2, 0.3, 0.5
	mat->transparency = materialSetting->getTransparency();																// 0.7
	res->addChild(mat);

	// Gouraud shading
	SoShapeHints* hints = new SoShapeHints;
	hints->creaseAngle = 3.14;
	res->addChild(hints);

	SoCoordinate3* coords = new SoCoordinate3();
	for (int c = 0; c < mesh->getNumOfVerts(); c++)
		coords->point.set1Value(c, mesh->getVertex(c).coords[0],
			mesh->getVertex(c).coords[1],
			mesh->getVertex(c).coords[2]);

	SoIndexedFaceSet* faceSet = new SoIndexedFaceSet();
	for (int c = 0; c < mesh->getNumOfTris(); c++)
	{
		faceSet->coordIndex.set1Value(c * 4, mesh->getTriangle(c).corners[0]);
		faceSet->coordIndex.set1Value(c * 4 + 1, mesh->getTriangle(c).corners[1]);
		faceSet->coordIndex.set1Value(c * 4 + 2, mesh->getTriangle(c).corners[2]);
		faceSet->coordIndex.set1Value(c * 4 + 3, -1);		
	}

	res->addChild(coords);
	res->addChild(faceSet);
}

void Painter::getShapeSep(double center[3], float radius, SoSeparator* res)
{
	// Paint all vertices with the same color
	SoMaterial* mat = new SoMaterial();
	mat->diffuseColor.setValue(0.5, 0.5, 0.5);
	mat->transparency = 0;
	res->addChild(mat);

	// Gouraud shading
	SoShapeHints* hints = new SoShapeHints;
	hints->creaseAngle = 3.14;
	res->addChild(hints);

	//transformation
	SoTransform* transform = new SoTransform();
	transform->translation.setValue((float*)center);
	res->addChild(transform);

	SoSphere* sphere = new SoSphere();
	sphere->radius = radius;
	res->addChild(sphere);
}

void Painter::getShapeSep(Mesh* mesh, vector<double*>& coveringNodes, vector<int*>& trisIds, vector<int>& numOfTrisPerNode, float radius, SoSeparator* res) {

	const int numOfColors = 12;
	double colors[numOfColors][3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 0}, {1, 0, 1}, {0, 1, 1}, 
						  {0.2, 0.5, 0.8}, {0.2, 0.8, 0.5}, {0.8, 0.5, 0.2}, {0.8, 0.2, 0.5}, {0.5, 0.2, 0.8}, {0.5, 0.8, 0.2}};
	vector<int> drawnTris;

	for (int i = 0; i < coveringNodes.size(); i++) {
		SoSeparator* segmentSep = new SoSeparator;

		// DRAW TRIANGLES
		SoSeparator* trisSep = new SoSeparator;
		// Paint all vertices with the same color
		SoMaterial* mat = new SoMaterial();
		mat->transparency = 0.5;
		mat->diffuseColor.setValue(colors[i % numOfColors][0], colors[i % numOfColors][1], colors[i % numOfColors][2]);
		trisSep->addChild(mat);

		// Add vertex info
		SoCoordinate3* coords = new SoCoordinate3();
		for (int c = 0; c < mesh->getNumOfVerts(); c++)
			coords->point.set1Value(c, mesh->getVertex(c).coords[0],
				mesh->getVertex(c).coords[1],
				mesh->getVertex(c).coords[2]);
		trisSep->addChild(coords);

		// Add triangle info
		SoIndexedFaceSet* faceSet = new SoIndexedFaceSet();
		for (int c = 0, t = 0; c < numOfTrisPerNode[i]; c++)
		{
			if (count(drawnTris.begin(), drawnTris.end(), trisIds[i][c]) == 0) {
				faceSet->coordIndex.set1Value(t * 4, mesh->getTriangle(trisIds[i][c]).corners[0]);
				faceSet->coordIndex.set1Value(t * 4 + 1, mesh->getTriangle(trisIds[i][c]).corners[1]);
				faceSet->coordIndex.set1Value(t * 4 + 2, mesh->getTriangle(trisIds[i][c]).corners[2]);
				faceSet->coordIndex.set1Value(t * 4 + 3, -1);
				drawnTris.push_back(trisIds[i][c]);
				t++;
			}
		}
		trisSep->addChild(faceSet);
		segmentSep->addChild(trisSep);

		// DRAW SPHERE
		SoSeparator* sphereSep = new SoSeparator;
		// Paint all vertices with the same color
		SoMaterial* mat2 = new SoMaterial();
		mat2->diffuseColor.setValue(colors[i % numOfColors][0], colors[i % numOfColors][1], colors[i % numOfColors][2]);
		sphereSep->addChild(mat2);

		//transformation
		SoTransform* transform = new SoTransform();
		float pcoords[3];
		for (int j = 0; j < 3; j++)
			pcoords[j] = coveringNodes[i][j];
		transform->translation.setValue(pcoords);
		sphereSep->addChild(transform);

		SoSphere* sphere = new SoSphere();
		sphere->radius = radius;
		sphereSep->addChild(sphere);
		segmentSep->addChild(sphereSep);

		res->addChild(segmentSep);
	}
}

void Painter::getShapeSepByEdges(Mesh* mesh, MaterialSetting* materialSetting, SoSeparator* res)
{
	// Paint all vertices with the same color
	SoMaterial* mat = new SoMaterial();
	mat->diffuseColor.setValue(materialSetting->getRed(), materialSetting->getGreen(), materialSetting->getBlue());		// 0.2, 0.3, 0.5
	mat->transparency = materialSetting->getTransparency();																// 0.7
	res->addChild(mat);

	// Gouraud shading
	SoShapeHints* hints = new SoShapeHints;
	hints->creaseAngle = 3.14;
	res->addChild(hints);

	SoCoordinate3* coords = new SoCoordinate3();
	for (int c = 0; c < mesh->getNumOfVerts(); c++)
		coords->point.set1Value(c, mesh->getVertex(c).coords[0],
			mesh->getVertex(c).coords[1],
			mesh->getVertex(c).coords[2]);

	SoDrawStyle* sty = new SoDrawStyle;
	sty->lineWidth = 1.0f;
	res->addChild(sty);
	SoIndexedLineSet* ilsSet = new SoIndexedLineSet;
	for (int c = 0; c < mesh->getNumOfEdges(); c++)
	{
		ilsSet->coordIndex.set1Value(c * 2, mesh->getEdge(c).endVerts[0]);
		ilsSet->coordIndex.set1Value(c * 2 + 1, mesh->getEdge(c).endVerts[1]);
		ilsSet->coordIndex.set1Value(2, -1);
	}

	res->addChild(coords);
	res->addChild(ilsSet);
}

void Painter::drawTriangulation(Mesh* mesh, SoSeparator* res) {

	SoSeparator* thickEdgeSep = new SoSeparator;
	SoMaterial* ma = new SoMaterial;
	//ma->diffuseColor.set1Value(0, 1.0, 0.4, 0.4);
	ma->diffuseColor.set1Value(0, 1.0, 0, 0);
	thickEdgeSep->addChild(ma);
	SoDrawStyle* sty = new SoDrawStyle;
	sty->lineWidth = 1.0f;
	thickEdgeSep->addChild(sty);

	SoIndexedLineSet* ils = new SoIndexedLineSet;
	SoCoordinate3* co = new SoCoordinate3;
	for (int se = 0; se < mesh->getNumOfEdges(); se++)
	{
		co->point.set1Value(2 * se, mesh->getVertex(mesh->getEdge(se).endVerts[0]).coords[0],
									mesh->getVertex(mesh->getEdge(se).endVerts[0]).coords[1], 
									mesh->getVertex(mesh->getEdge(se).endVerts[0]).coords[2]);
		co->point.set1Value(2 * se + 1, mesh->getVertex(mesh->getEdge(se).endVerts[1]).coords[0],
										mesh->getVertex(mesh->getEdge(se).endVerts[1]).coords[1], 
										mesh->getVertex(mesh->getEdge(se).endVerts[1]).coords[2]);
	}

	for (int ci = 0; ci < mesh->getNumOfEdges(); ci++)
	{
		ils->coordIndex.set1Value(3 * ci, 2 * ci);
		ils->coordIndex.set1Value(3 * ci + 1, 2 * ci + 1);
		ils->coordIndex.set1Value(3 * ci + 2, -1);
	}

	thickEdgeSep->addChild(co);
	thickEdgeSep->addChild(ils);
	res->addChild(thickEdgeSep);

}

void Painter::drawColorfulTriangulation(Mesh* mesh, SoSeparator* res) {

	SoSeparator* thickEdgeSep = new SoSeparator;
	SoMaterial* ma = new SoMaterial;
	thickEdgeSep->addChild(ma);

	SoCoordinate3* co = new SoCoordinate3();
	for (int c = 0; c < mesh->getNumOfVerts(); c++) {
		co->point.set1Value(c, mesh->getVertex(c).coords[0], mesh->getVertex(c).coords[1], mesh->getVertex(c).coords[2]);
		ma->diffuseColor.set1Value(c, mesh->getVertex(c).color[0], mesh->getVertex(c).color[1], mesh->getVertex(c).color[2]);
	}

	SoMaterialBinding* materialBinding = new SoMaterialBinding; //for 2+ diffuse color usage on the same mesh
	materialBinding->value = SoMaterialBinding::PER_VERTEX_INDEXED;
	thickEdgeSep->addChild(materialBinding);
	thickEdgeSep->addChild(co);

	SoDrawStyle* sty = new SoDrawStyle;
	sty->lineWidth = 1.0f;
	thickEdgeSep->addChild(sty);
	SoIndexedLineSet* ils = new SoIndexedLineSet;

	for (int ci = 0; ci < mesh->getNumOfEdges(); ci++)
	{
		ils->coordIndex.set1Value(3 * ci, mesh->getEdge(ci).endVerts[0]);
		ils->coordIndex.set1Value(3 * ci + 1, mesh->getEdge(ci).endVerts[1]);
		ils->coordIndex.set1Value(3 * ci + 2, -1);
	}

	thickEdgeSep->addChild(ils);
	res->addChild(thickEdgeSep);

}

void Painter::drawSingleTriangle(Mesh* mesh, Triangle* triangle, SoSeparator* res) {

	SoSeparator* thickEdgeSep = new SoSeparator;
	SoMaterial* ma = new SoMaterial;
	ma->diffuseColor.set1Value(0, 0.0f, 0.0f, 1.0f);	// blue
	thickEdgeSep->addChild(ma);
	SoDrawStyle* sty = new SoDrawStyle;
	sty->lineWidth = 1.0f;
	thickEdgeSep->addChild(sty);

	SoIndexedLineSet* ils = new SoIndexedLineSet;
	SoCoordinate3* co = new SoCoordinate3;

	for (int se = 0; se < 3; se++) {
		Edge edge = mesh->getEdge(triangle->edgeList[se]);
		co->point.set1Value(2 * se, mesh->getVertex(edge.endVerts[0]).coords[0],
			mesh->getVertex(edge.endVerts[0]).coords[1],
			mesh->getVertex(edge.endVerts[0]).coords[2]);
		co->point.set1Value(2 * se + 1, mesh->getVertex(edge.endVerts[1]).coords[0],
			mesh->getVertex(edge.endVerts[1]).coords[1],
			mesh->getVertex(edge.endVerts[1]).coords[2]);
	}

	for (int se = 0; se < 3; se++) {
		ils->coordIndex.set1Value(3 * se, 2 * se);
		ils->coordIndex.set1Value(3 * se + 1, 2 * se + 1);
		ils->coordIndex.set1Value(3 * se + 2, -1);
	}

	thickEdgeSep->addChild(co);
	thickEdgeSep->addChild(ils);
	res->addChild(thickEdgeSep);

}

void Painter::drawNormal(Mesh* mesh, int triangle_id, SoSeparator* res) {

	Triangle triangle = mesh->getTriangle(triangle_id);
	SoMaterial* mat = new SoMaterial();
	mat->diffuseColor.setValue(0.5, 0.5, 0.5);
	mat->transparency = 0;
	res->addChild(mat);

	// compute centerx3
	double centerx3[3] = { 0, 0, 0 };
	for (int i = 0; i < 3; i++) {
		int vertex_id = triangle.corners[i];
		for (int j = 0; j < 3; j++)
			centerx3[j] += mesh->getVertex(vertex_id).coords[j];
	}

	// compute center and the other end of the normal
	SoCoordinate3* coords = new SoCoordinate3();
	coords->point.set1Value(0, centerx3[0] / 3, centerx3[1] / 3, centerx3[2] / 3);
	coords->point.set1Value(1, centerx3[0] / 3 + (triangle.normal[0] / 100),
		centerx3[1] / 3 + (triangle.normal[1] / 100), centerx3[2] / 3 + (triangle.normal[2] / 100));

	// Draw normal
	SoDrawStyle* sty = new SoDrawStyle;
	sty->lineWidth = 1.0f;
	res->addChild(sty);
	SoIndexedLineSet* ils = new SoIndexedLineSet;
	ils->coordIndex.set1Value(0, 0);
	ils->coordIndex.set1Value(1, 1);
	ils->coordIndex.set1Value(2, -1);

	res->addChild(coords);
	res->addChild(ils);
}

void Painter::drawLine(double* p1, double* p2, SoSeparator* res) {

	SoMaterial* mat = new SoMaterial();
	mat->diffuseColor.setValue(0, 0, 1.0);
	mat->transparency = 0;
	res->addChild(mat);

	SoCoordinate3* coords = new SoCoordinate3();
	coords->point.set1Value(0, p1[0], p1[1], p1[2]);
	coords->point.set1Value(1, p2[0], p2[1], p2[2]);

	// Draw normal
	SoDrawStyle* sty = new SoDrawStyle;
	sty->lineWidth = 5.0f;
	res->addChild(sty);
	SoIndexedLineSet* ils = new SoIndexedLineSet;
	ils->coordIndex.set1Value(0, 0);
	ils->coordIndex.set1Value(1, 1);
	ils->coordIndex.set1Value(2, -1);

	res->addChild(coords);
	res->addChild(ils);
}


