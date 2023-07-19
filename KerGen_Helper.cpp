#include "KerGen.h"
#include "BaseGeoOpUtils.h"

bool KerGen::isLineDirectionValid(double* startPoint, double* lineDirection) {

	double testPoint[3];
	for (int k = 0; k < 3; k++)
		testPoint[k] = startPoint[k] + lineDirection[k] * EPSILON;

	double* scalars = findValueVectorSatisfiedByPoint(testPoint, halfSpaceSet);

	bool isValid = true;
	for (int k = 0; k < halfSpaceSet.size(); k++)
		if (scalars[k] > EPSILON) {
			isValid = false;
			break;
		}
	delete[] scalars;

	return isValid;
}

bool KerGen::isLineValid(double* startPoint, double* lineDirection) {

	double testPoint[3];
	for (int k = 0; k < 3; k++)
		testPoint[k] = startPoint[k] + lineDirection[k] * EPSILON;

	double* scalars = findValueVectorSatisfiedByPoint(testPoint, halfSpaceSet);

	bool isValid = true;
	for (int k = 0; k < halfSpaceSet.size(); k++)
		if (scalars[k] > 3 * EPSILON) {
			isValid = false;
			break;
		}
	delete[] scalars;

	if (!isValid) {	// revert the line direction
		for (int k = 0; k < 3; k++) {
			lineDirection[k] *= -1;
			testPoint[k] = startPoint[k] + lineDirection[k] * EPSILON;
		}

		double* scalars = findValueVectorSatisfiedByPoint(testPoint, halfSpaceSet);

		isValid = true;
		for (int k = 0; k < halfSpaceSet.size(); k++)
			if (scalars[k] > 3 * EPSILON) {
				isValid = false;
				break;
			}
		delete[] scalars;
	}

	return isValid;
}

bool KerGen::shiftOnAnotherPlane(double* point, double* rootpoint, int& shiftId) {

	bool shouldShift = false;
	double* scalars = findValueVectorSatisfiedByPoint(point, halfSpaceSet);
	double theInnerMostDistance = EPSILON;
	int theInnerMostId = -1;
	for (int i = 0; i < halfSpaceSet.size(); i++)
		if (scalars[i] > theInnerMostDistance) {

			shouldShift = true;
			if (isPointOnPlane(rootpoint, i)) {
				theInnerMostDistance = scalars[i];
				theInnerMostId = i;
			}

		}
	delete[] scalars;

	shiftId = theInnerMostId;
	return shouldShift;

}

bool KerGen::shiftOnAnotherPlane(double* point, double* rootpoint, int& shiftId, int baseId) {

	bool shouldShift = false;
	double* scalars = findValueVectorSatisfiedByPoint(point, halfSpaceSet);
	double theInnerMostDistance = -1;
	int theInnerMostId = -1;
	for (int i = 0; i < halfSpaceSet.size(); i++) {
		if (scalars[i] > EPSILON) {
			shouldShift = true;
			if (isPointOnPlane(rootpoint, i)) {
				double lineDir[3];
				crossProduct(halfSpaceSet[baseId].ABCD, halfSpaceSet[i].ABCD, lineDir);
				if (isZeroVector(lineDir))
					continue;
				normalize(lineDir);
				double rayDir[3];
				crossProduct(halfSpaceSet[baseId].ABCD, lineDir, rayDir);
				normalize(rayDir);

				Line ray(point, rayDir);
				double t = findLinePlaneIntersection(ray, halfSpaceSet[i]);

				if (abs(t) > theInnerMostDistance) {
					theInnerMostDistance = abs(t);
					theInnerMostId = i;
				}
			}
		}
	}
	delete[] scalars;

	shiftId = theInnerMostId;
	return shouldShift;

}

bool KerGen::isPointOnPlane(double* point, int planeId) {

	double distance = halfSpaceSet[planeId].ABCD[3];
	for (int d = 0; d < 3; d++)
		distance += halfSpaceSet[planeId].ABCD[d] * point[d];

	if (abs(distance) < EPSILON)
		return true;
	return false;

}

bool KerGen::isPointInFrontOfThePlane(double* point, int planeId) {

	double distance = halfSpaceSet[planeId].ABCD[3];
	for (int d = 0; d < 3; d++)
		distance += halfSpaceSet[planeId].ABCD[d] * point[d];

	if (distance < EPSILON)
		return true;
	return false;

}

bool KerGen::areTheSamePlanes(int plane1Id, int plane2Id) {

	if (isTripleSame(halfSpaceSet[plane1Id].ABCD, halfSpaceSet[plane2Id].ABCD)) {
		double* point = halfSpaceSet[plane1Id].point1;
		if (isPointOnPlane(point, plane2Id))
			return true;
	}

	return false;

}

void KerGen::findInitialVertexAndLines(double* point) {

	double distanceToFirst;
	vector<int> firstClosestIds = findFirstPlane(point, distanceToFirst);

	for (int i = 0; i < firstClosestIds.size(); i++) {
		int plane1Id = firstClosestIds[i];

		for (int d = 0; d < 3; d++)
			point[d] = point[d] + halfSpaceSet[plane1Id].ABCD[i] * distanceToFirst;

		vector<double> pointerToLineDirections, intersectionLineDirections;
		double distanceToSecond;
		vector<int> secondClosestIds = findSecondPlane(point, plane1Id, pointerToLineDirections, intersectionLineDirections, distanceToSecond);

		for (int j = 0; j < secondClosestIds.size(); j++) {
			int plane2Id = secondClosestIds[j];

			double pointerToLineDir[3], firstLineDir[3];
			for (int d = 0; d < 3; d++) {
				pointerToLineDir[d] = pointerToLineDirections[j * 3 + d];
				firstLineDir[d] = intersectionLineDirections[j * 3 + d];
			}

			for (int d = 0; d < 3; d++)
				point[d] = point[d] + pointerToLineDir[d] * distanceToSecond;

			double vertex[3];
			vector<int> thirdClosestIds = findThirdPlane(point, firstLineDir, plane1Id, plane2Id, vertex);

			if (identifyInitialPlanes(vertex))
				return;

		}
	}

}

bool KerGen::identifyInitialPlanes(double* vertex) {

	vector<int> planeIds = vertexParentIds[0];
	//cout << planeIds.size() << endl;

	for (int i = 0; i < planeIds.size(); i++) {
		int planeId1 = planeIds[i];
		bool plane1_valid = false;

		for (int j = i; j < planeIds.size(); j++) {
			int planeId2 = planeIds[j];
			bool plane2_valid = false;

			if (isTripleSame(halfSpaceSet[planeId1].ABCD, halfSpaceSet[planeId2].ABCD))
				continue;

			double lineDir1[3];
			crossProduct(halfSpaceSet[planeId1].ABCD, halfSpaceSet[planeId2].ABCD, lineDir1);
			normalize(lineDir1);

			if (!plane1_valid) {
				if (isLineValid(vertex, lineDir1)) {
					double perpDir1[3];
					crossProduct(halfSpaceSet[planeId1].ABCD, lineDir1, perpDir1);
					normalize(perpDir1);
					double point[3];
					for (int d = 0; d < 3; d++)
						point[d] = vertex[d] + EPSILON * lineDir1[d];
					if (isLineValid(point, perpDir1))
						plane1_valid = true;
					else
						continue;
				}
				else
					continue;
			}

			double perpDir2[3];
			crossProduct(halfSpaceSet[planeId2].ABCD, lineDir1, perpDir2);
			normalize(perpDir2);
			double point[3];
			for (int d = 0; d < 3; d++)
				point[d] = vertex[d] + EPSILON * lineDir1[d];
			if (isLineValid(point, perpDir2))
				plane2_valid = true;
			else
				continue;

			for (int k = j; k < planeIds.size(); k++) {
				int planeId3 = planeIds[k];
				bool plane3_valid = false;

				if (isTripleSame(halfSpaceSet[planeId1].ABCD, halfSpaceSet[planeId3].ABCD))
					continue;
				if (isTripleSame(halfSpaceSet[planeId2].ABCD, halfSpaceSet[planeId3].ABCD))
					continue;

				double lineDir2[3];
				crossProduct(halfSpaceSet[planeId1].ABCD, halfSpaceSet[planeId3].ABCD, lineDir2);
				normalize(lineDir2);

				if (isLineValid(vertex, lineDir2)) {
					double perpDir3[3];
					crossProduct(halfSpaceSet[planeId1].ABCD, lineDir2, perpDir3);
					normalize(perpDir3);
					double point[3];
					for (int d = 0; d < 3; d++)
						point[d] = vertex[d] + EPSILON * lineDir2[d];
					if (isLineValid(point, perpDir3))
						plane3_valid = true;
				}

				double lineDir3[3];
				crossProduct(halfSpaceSet[planeId2].ABCD, halfSpaceSet[planeId3].ABCD, lineDir3);
				normalize(lineDir3);

				if (isLineValid(vertex, lineDir3)) {
					double perpDir3[3];
					crossProduct(halfSpaceSet[planeId2].ABCD, lineDir3, perpDir3);
					normalize(perpDir3);
					double point[3];
					for (int d = 0; d < 3; d++)
						point[d] = vertex[d] + EPSILON * lineDir3[d];
					if (isLineValid(point, perpDir3))
						plane3_valid = true;
					else
						continue;
				}
				else
					continue;

				// if I come here, I successfully found 3 planes
				recordLineInfo(planeId1, planeId2, lineDir1, vertex);
				recordLineInfo(planeId1, planeId3, lineDir2, vertex);
				recordLineInfo(planeId2, planeId3, lineDir3, vertex);

				return true;

			}

		}

	}

	return false;
}

vector<int> KerGen::findFirstPlane(double* point, double& distanceToFirst) {

	double* scalars = findValueVectorSatisfiedByPoint(point, halfSpaceSet);

	distanceToFirst = numeric_limits<double>::infinity();
	vector<int> theClosestIds;
	for (int i = 0; i < halfSpaceSet.size(); i++) {
		if (abs(scalars[i]) <= distanceToFirst + BIG_EPSILON) {
			if (abs(scalars[i]) < distanceToFirst - BIG_EPSILON)
				theClosestIds.clear();
			theClosestIds.push_back(i);
			distanceToFirst = abs(scalars[i]);
		}
	}
	delete[] scalars;

	cout << theClosestIds.size() << endl;
	return theClosestIds;

}

vector<int> KerGen::findSecondPlane(double* point, int id, vector<double>& pointerToLineDirections, vector<double>& intersectionLineDirections, double& distanceToSecond) {

	distanceToSecond = numeric_limits<double>::infinity();
	int theClosestId = -1;
	vector<int> theClosestIds;
	for (int s = 0, i = 0, bound = id; s < 2; s++) {
		for (; i < bound; i++) {
			if (isTripleSame(halfSpaceSet[id].ABCD, halfSpaceSet[i].ABCD))
				continue;
			double lineDirection[3], rayDirection[3];
			crossProduct(halfSpaceSet[i].ABCD, halfSpaceSet[id].ABCD, lineDirection);
			normalize(lineDirection);
			crossProduct(halfSpaceSet[id].ABCD, lineDirection, rayDirection);
			normalize(rayDirection);
			Line ray(point, rayDirection);
			double t = findLinePlaneIntersection(ray, halfSpaceSet[i]);
			if (abs(t) <= distanceToSecond + BIG_EPSILON) {
				if (abs(t) < distanceToSecond - BIG_EPSILON) {
					theClosestIds.clear();
					pointerToLineDirections.clear();
					intersectionLineDirections.clear();
				}
				theClosestIds.push_back(i);
				distanceToSecond = abs(t);
				for (int k = 0; k < 3; k++) {
					pointerToLineDirections.push_back(rayDirection[k]);
					intersectionLineDirections.push_back(lineDirection[k]);
				}
			}
		}
		i = id + 1;
		bound = halfSpaceSet.size();
	}

	return theClosestIds;

}

vector<int> KerGen::findThirdPlane(double* startpoint, double* lineDirection, int lineParent1Id, int lineParent2Id, double* newpoint) {

	double theClosestDistance = numeric_limits<double>::infinity();
	vector<int> theClosestIds;
	vector<int> theSameOfLineParents;

	Line line(startpoint, lineDirection);

	for (int i = 0; i < halfSpaceSet.size(); i++) {

		if (areTheSamePlanes(lineParent1Id, i)) {
			theSameOfLineParents.push_back(i);
			continue;
		}
		if (areTheSamePlanes(lineParent2Id, i)) {
			theSameOfLineParents.push_back(i);
			continue;
		}

		double t = findLinePlaneIntersection(line, halfSpaceSet[i]);

		if (t >= -BIG_EPSILON && t <= theClosestDistance + BIG_EPSILON) {
			if (t < theClosestDistance - BIG_EPSILON)
				theClosestIds.clear();
			theClosestIds.push_back(i);
			theClosestDistance = t;
		}

	}

	for (int k = 0; k < 3; k++)
		newpoint[k] = startpoint[k] + lineDirection[k] * theClosestDistance;

	recordVertexInfo(newpoint, theClosestIds, theSameOfLineParents);

	return theClosestIds;

}

bool KerGen::isPreviousLine(int hs1_id, int hs2_id) {

	bool walkedEdge = false;
	for (int j = 0; j < edgePartnerIds[hs1_id].size(); j++)
		if (edgePartnerIds[hs1_id][j] == hs2_id) {
			walkedEdge = true;
			break;
		}
	return walkedEdge;
}

void KerGen::recordLineInfo(int hs1_id, int hs2_id, double* lineDir, double* point) {

	edgePartners[hs1_id].push(EdgePartnerTriple(hs2_id, lineDir, point, kernel.getNumOfVerts() - 1));
	edgePartnerIds[hs1_id].push_back(hs2_id);
	edgePartnerIds[hs2_id].push_back(hs1_id);
	if (isKernelFace[hs1_id] != ProcessColor::GREEN)
		isKernelFace[hs1_id] = ProcessColor::WHITE;
	if (isKernelFace[hs2_id] != ProcessColor::GREEN)
		isKernelFace[hs2_id] = ProcessColor::WHITE;

}

int KerGen::recordVertexInfo(double* vertex, vector<int> thirdParentIds, vector<int> lineParentIds) {

	int vertexId = -1;
	for (int v = 0; v < kernel.getNumOfVerts(); v++)
		if (isTripleSame(kernel.getVertex(v).coords, vertex)) {
			vertexId = v;
			break;
		}

	if (vertexId == -1) {
		vertexId = kernel.getNumOfVerts();
		kernel.addVertex(vertex[0], vertex[1], vertex[2]);
		vector<int> parentIdsForNewVertex;
		vertexParentIds.push_back(parentIdsForNewVertex);

		for (int i = 0; i < thirdParentIds.size(); i++)
			vertexParentIds[vertexId].push_back(thirdParentIds[i]);
		for (int i = 0; i < lineParentIds.size(); i++)
			vertexParentIds[vertexId].push_back(lineParentIds[i]);
	}

	return vertexId;
}

vector<int> KerGen::findClosestPlaneOnLine(int startVertexId, double* lineDirection, int lineParent1Id, int lineParent2Id, double* vertex) {

	double theClosestDistance = numeric_limits<double>::infinity();
	vector<int> theClosestIds;
	double* startpoint = kernel.getVertex(startVertexId).coords;
	vector<int> theSameOfLineParents;

	while (true) {
		Line line(startpoint, lineDirection);

		for (int i = 0; i < halfSpaceSet.size(); i++) {


			if (areTheSamePlanes(lineParent1Id, i)) {
				theSameOfLineParents.push_back(i);
				continue;
			}
			if (areTheSamePlanes(lineParent2Id, i)) {
				theSameOfLineParents.push_back(i);
				continue;
			}

			
			double t = findLinePlaneIntersection(line, halfSpaceSet[i]);
			if (abs(t) < EPSILON)
				continue;
			
			/*
			bool same_of_parent = false;
			for (int p = 0; p < vertexParentIds[startVertexId].size(); p++)
				if (i == vertexParentIds[startVertexId][p]) {
					same_of_parent = true;
					break;
				}
			if (same_of_parent)
				continue;

			double t = findLinePlaneIntersection(line, halfSpaceSet[i]);
			*/
			if (t >= -EPSILON && t <= theClosestDistance + EPSILON) {
				if (t < theClosestDistance - EPSILON)
					theClosestIds.clear();
				theClosestIds.push_back(i);
				theClosestDistance = t;
			}
		}

		for (int d = 0; d < 3; d++)
			vertex[d] = startpoint[d] + lineDirection[d] * theClosestDistance;

		int shiftId;
		shiftOnAnotherPlane(vertex, startpoint, shiftId);
		//shiftOnAnotherPlane(vertex, startpoint, shiftId, lineParent1Id);
		if (shiftId >= 0) {
			//cout << "!!!!!!!!!!!!!!!!!\n";
			lineParent2Id = shiftId;
			crossProduct(halfSpaceSet[lineParent1Id].ABCD, halfSpaceSet[lineParent2Id].ABCD, lineDirection);
			normalize(lineDirection);
			isLineValid(startpoint, lineDirection);

			//theClosestDistance = numeric_limits<double>::infinity();
			//theClosestIds.clear();
			continue;	// repeat the steps with the new line
		}

		break;
	}

	int vertexId = recordVertexInfo(vertex, theClosestIds, theSameOfLineParents);

	kernel.addEdge(startVertexId, vertexId);

	return theClosestIds;

}

vector<int> KerGen::findClosestPlaneOnLine2(int startVertexId, double* lineDirection, int lineParent1Id, int lineParent2Id, double* vertex) {

	double theClosestDistance = numeric_limits<double>::infinity();
	vector<int> theClosestIds;
	double* startpoint = kernel.getVertex(startVertexId).coords;
	vector<int> theSameOfLineParents;

	while (true) {
		Line line(startpoint, lineDirection);

		for (int i = 0; i < halfSpaceSet.size(); i++) {

			if (areTheSamePlanes(lineParent1Id, i)) {
				theSameOfLineParents.push_back(i);
				continue;
			}
			if (areTheSamePlanes(lineParent2Id, i)) {
				theSameOfLineParents.push_back(i);
				continue;
			}

			bool previous_parent = false;
			for (int j = 0; j < vertexParentIds[startVertexId].size(); j++)
				if (i == vertexParentIds[startVertexId][j]) {
					previous_parent = true;
					break;
				}
				else if (areTheSamePlanes(i, vertexParentIds[startVertexId][j])) {
					previous_parent = true;
					vertexParentIds[startVertexId].push_back(i);
					break;
				}
			if (previous_parent)
				continue;

			double t = findLinePlaneIntersection(line, halfSpaceSet[i]);

			if (t >= 0 && t <= theClosestDistance) {
				if (t < theClosestDistance - EPSILON)
					theClosestIds.clear();
				theClosestIds.push_back(i);
				theClosestDistance = t;
			}
		}

		for (int d = 0; d < 3; d++)
			vertex[d] = startpoint[d] + lineDirection[d] * theClosestDistance;

		int shiftId;
		shiftOnAnotherPlane(vertex, startpoint, shiftId, lineParent1Id);
		if (shiftId >= 0) {
			theClosestDistance = numeric_limits<double>::infinity();
			theClosestIds.clear();
			lineParent2Id = shiftId;
			crossProduct(halfSpaceSet[lineParent1Id].ABCD, halfSpaceSet[lineParent2Id].ABCD, lineDirection);
			normalize(lineDirection);
			isLineValid(startpoint, lineDirection);
			continue;	// repeat the steps with the new line
		}

		break;
	}

	int vertexId = recordVertexInfo(vertex, theClosestIds, theSameOfLineParents);

	kernel.addEdge(startVertexId, vertexId);

	return theClosestIds;

}

