//
// Created by Robert on 10/28/2019.
//

#ifndef MINACALC_SMLOADER_H
#define MINACALC_SMLOADER_H

#include <vector>
#include <fstream>
#include "NoteDataStructures.h"

using std::ifstream;

typedef std::vector<NoteInfo> SMNotes;

SMNotes load_from_file(ifstream& sm_file);

#endif //MINACALC_SMLOADER_H
