#ifndef _MEDIA_READ_INTERFACE_FILE_
#define _MEDIA_READ_INTERFACE_FILE_
#define MAX_BUF_LEN 4096

#include <iostream>
#include <vector>
#include "media_utility.hpp"

FILE *gfopen(const char *filename, const char *mode);

void read_interface_file(
    const char *interface_file,
    inter_t **interfaces,
    int &md_type);


void read_grid_file(
    const char *grid_file,
    float Xmin, float Xmax,
    int &NL,
    std::vector<int> &NGz, // how many z-grid in each layer
    inter_t *interfaces);

// check whether the elevation[ng[i]-1] == elevation[ng[i]] 
int checkGridData(int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces, 
    const char *grid_file); 

#endif