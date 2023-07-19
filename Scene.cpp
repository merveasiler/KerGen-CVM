// @author Merve Asiler

#include "Scene.h"
#include <iostream>

SoPath* pickFilterCB(void*, const SoPickedPoint* pick)
{
	// See which child of selection got picked
	SoPath* p = pick->getPath();
	int i;
	for (i = 0; i < p->getLength() - 1; i++) {
		SoNode* n = p->getNode(i);
		if (n->isOfType(SoSelection::getClassTypeId()))
			break;
	}

	// Copy 2 nodes from the path:
	// selection and the picked child
	return p->copy(i, 2);
}


Scene::Scene() {

	window = SoWin::init("Kernel Computation");
	viewer = new SoWinExaminerViewer(window);
	root = new SoSeparator;
	root->ref();

}

void Scene::attachToRoot(SoSeparator* res) {
	root->addChild(res);
}

void Scene::attachToRoot(SoSwitch* res) {
	root->addChild(res);
}

void Scene::configureSelection(SoSeparator* res) {

	selection = new SoSelection;
	selection->addChild(res);
	selection->setPickFilterCallback(pickFilterCB);

}

void Scene::configureViewer() {

	viewer->setBackgroundColor(SbColor(1, 1, 1));
	viewer->setSize(SbVec2s(640, 480));
	viewer->setSceneGraph(root);
	viewer->viewAll();
	viewer->show();
	//viewer->setSceneGraph(selection);
}

SoOrthographicCamera* Scene::configureOrthographicCamera(float cameraHeight) {

	SoOrthographicCamera* camera = new SoOrthographicCamera;
	camera->position.setValue(0, 0, 25.0);
	camera->pointAt(SbVec3f(0, 0, 1.0), SbVec3f(0, 1.0, 0));	// orientation
	camera->nearDistance.setValue(10.0);
	camera->farDistance.setValue(40.0);
	camera->focalDistance.setValue(25.0);
	camera->aspectRatio.setValue(1.0);							// SO_ASPECT_SQUARE
	camera->height = cameraHeight;

	return camera;

}

void Scene::play() {

	SoWin::show(window);
	SoWin::mainLoop();
}

Scene::~Scene() {
	root->unref();
	//selection->removeAllChildren();
	delete viewer;
}

void Scene::makeScene(SoSeparator* res) {

	attachToRoot(res);
	//configureSelection(res);
	configureViewer();
	play();
}

void Scene::makeScene(vector< SoSeparator* > resSet) {

	for (auto& res : resSet)
		attachToRoot(res);
	configureViewer();
	play();
}

void Scene::makeScene(vector< SoSwitch* > resSet) {

	for (auto& res : resSet)
		attachToRoot(res);
	configureViewer();
	play();
}

void Scene::makeMultipleScene(vector<SoSeparator*> resSets, SoSeparator* cornerResSet, float sceneHeight) {

	SoOrthographicCamera* camera1 = configureOrthographicCamera(sceneHeight);
	SoOrthographicCamera* camera2 = camera1;// configureOrthographicCamera(sceneHeight);

	if (cornerResSet != NULL) {
		root->addChild(camera1);
		root->addChild(cornerResSet);
		configureViewer();
	}

	root = new SoSeparator;
	root->addChild(camera2);
	for (int i=0; i < resSets.size(); i++)
		root->addChild(resSets[i]);
	root->ref();
	viewer = new SoWinExaminerViewer(window);
	configureViewer();

	play();

}



