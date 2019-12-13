//
// Created by Robert on 10/28/2019.
//

#ifndef MINACALC_SMLOADER_H
#define MINACALC_SMLOADER_H

#include <vector>
#include <fstream>
#include "NoteDataStructures.h"

typedef std::vector<ChartInfo> SMNotes;

struct BPM {
    float beat;
    float bpm;
};

typedef std::vector<BPM> BPMs;

SMNotes load_from_file(std::ifstream& sm_file);

#endif //MINACALC_SMLOADER_H
