#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>

using namespace std;
using namespace boost::filesystem;

void convertFromSTLToOFF(string inputFolder, string outputFolder);
