//
// Created by Robert on 10/28/2019.
//

#include <sstream>
#include <iostream>
#include "smloader.h"

using std::stringstream;
using std::string;
SMNotes parse_main_block(stringstream&);
BPMs parse_bpms_block(stringstream&);

SMNotes load_from_file(ifstream& file) {
    stringstream sm_buffer;
    sm_buffer << file.rdbuf();
    string sm_text = sm_buffer.str();
    SMNotes raw_block = std::vector<NoteInfo>();
    BPMs bpms = std::vector<BPM>();
    while (!sm_text.empty()) {
        size_t cool = sm_text.find('#');
        if (cool == string::npos)
            break;
        sm_text = sm_text.substr(cool+1);
        if (sm_text.substr(0,5) == "NOTES") {
            for (int i = 0; i < 6; i++) {
                cool = sm_text.find(':');
                if (cool == string::npos)
                    break;
                sm_text = sm_text.substr(cool+1);
            }
            stringstream notes_block;
            cool = sm_text.find(';');
            notes_block << sm_text.substr(sm_text.find('\n')+1, cool);
            raw_block = parse_main_block(notes_block);
            sm_text = sm_text.substr(cool+1);
        } else if (sm_text.substr(0,4) == "BPMS") {
            cool = sm_text.find(':');
            sm_text = sm_text.substr(cool+1);
            cool = sm_text.find(';');
            stringstream bpms_block;
            bpms_block << sm_text.substr(0, cool);
            bpms = parse_bpms_block(bpms_block);
            sm_text = sm_text.substr(cool+1);
        }
    }
    return raw_block;
}

SMNotes parse_main_block(stringstream& sm_text) {
    SMNotes output;
    int notes;
    int stupidname;
    float size = 0.f;
    float measurenumber = 0.f;
    std::vector<int> measure;
    for (std::string line; std::getline(sm_text, line); )
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
            if (stupidname == 8) {
                break;
            }
            stupidname *= 2;
        }
        size += 1.f;
        measure.push_back(notes);
    }
    float inside = 0.f;
    for(int note_row : measure) {
        if (note_row != 0) {
            output.push_back(NoteInfo{note_row, measurenumber + inside / size});
        }
        inside += 1.f;
    }
    measure.clear();

    return output;
}

BPMs parse_bpms_block(stringstream& bpms_block) {
    float next_time;
    float next_bpm;
    BPMs bpm_list = std::vector<BPM>();
    while (bpms_block.rdbuf()->in_avail() && (bpms_block >> next_time) && (bpms_block.get() == '=') && (bpms_block >> next_bpm)) {
        bpm_list.push_back(BPM{next_time, next_bpm});
        bpms_block.get();
        bpms_block.get();
    }
    return bpm_list;
}