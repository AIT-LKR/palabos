#include <memory>
#include <string>
#include <sstream>
#include <vector>
#include <chrono>
#include <random>
#include <sys/stat.h>

#ifndef TOOLS_FILES_HH
#define TOOLS_FILES_HH

using namespace std;


bool doesFileExist(std::string fileName) {
    struct stat buffer;
    return (stat (fileName.c_str(), &buffer) == 0);
}

bool checkCreateDir(std::string dirName) {
    bool created;
    // https://stackoverflow.com/questions/7430248/creating-a-new-directory-in-c
    if (doesFileExist(dirName)) {
        created = false;
    } else {
        mkdir(dirName.c_str(), 0755);
        created = true;
    }
    return created;
}


#endif
