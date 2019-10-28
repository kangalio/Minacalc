//
// Created by Robert on 10/28/2019.
//

#include <sstream>
#include "smloader.h"

using std::stringstream;
using std::string;

SMNotes load_from_file(ifstream& file) {
    stringstream sm_buffer;
    SMNotes output;
    sm_buffer << file.rdbuf();
    int notes;
    int stupidname;
    float size = 0.f;
    float measurenumber = 0.f;
    std::vector<int> measure;
    for (std::string line; std::getline(sm_buffer, line); )
    {
        notes = 0;
        stupidname = 1;
        if (line[0] == ',') {
            float inside = 0.f;
            for(int note_row : measure) {
                if (note_row != 0) {
                    output.push_back(NoteInfo {note_row, measurenumber + inside / size});
                }
                inside += 1.f;
            }
            measurenumber += 1.f;
            measure.clear();
            size = 0.f;
            continue;
        }
        for(char & it : line) {
            if (it == '1' || it == '2') {
                notes += stupidname;
            }
            stupidname *= 2;
        }
        size += 1.f;
        measure.push_back(notes);
    }

    return output;
}
