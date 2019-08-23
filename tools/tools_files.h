#include <memory>
#include <string>
#include <sstream>
#include <vector>
#include <chrono>
#include <random>
#include <sys/stat.h>

#ifndef TOOLS_FILES_H
#define TOOLS_FILES_H

using namespace std;


bool doesFileExist(std::string fileName);

bool checkCreateDir(std::string dirName);


#endif
