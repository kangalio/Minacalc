//
// Created by Robert on 10/28/2019.
//

#include <sstream>
#include <iostream>
#include <cctype>
#include "smloader.h"

using std::vector;
using std::stringstream;
using std::string;
vector<NoteInfo> parse_main_block(stringstream&);
BPMs parse_bpms_block(stringstream&);

SMNotes load_from_file(std::ifstream& file) {
    stringstream sm_buffer;
    sm_buffer << file.rdbuf();
    string sm_text = sm_buffer.str();
    SMNotes raw_block;
    BPMs bpms;
    while (!sm_text.empty()) {
        size_t next_tag_position = sm_text.find('#');
        if (next_tag_position == string::npos)
            break;
        sm_text = sm_text.substr(next_tag_position + 1);
        if (sm_text.substr(0,5) == "NOTES") {
            for (int i = 0; i < 3; i++) {
                next_tag_position = sm_text.find(':');
                if (next_tag_position == string::npos)
                    break;
                sm_text = sm_text.substr(next_tag_position + 1);
            }
            std::string difficulty_name = sm_text.substr(0,sm_text.find(':'));
            for (int i = 0; i < 3; i++) {
                next_tag_position = sm_text.find(':');
                if (next_tag_position == string::npos)
                    break;
                sm_text = sm_text.substr(next_tag_position + 1);
            }
            stringstream notes_block;
            next_tag_position = sm_text.find(';');
            notes_block << sm_text.substr(sm_text.find('\n')+1, next_tag_position-1);
            raw_block.push_back(ChartInfo {difficulty_name, parse_main_block(notes_block)});
            sm_text = sm_text.substr(next_tag_position + 1);
        } else if (sm_text.substr(0,4) == "BPMS") {
            next_tag_position = sm_text.find(':');
            sm_text = sm_text.substr(next_tag_position + 1);
            next_tag_position = sm_text.find(';');
            stringstream bpms_block;
            bpms_block << sm_text.substr(0, next_tag_position);
            bpms = parse_bpms_block(bpms_block);
            sm_text = sm_text.substr(next_tag_position + 1);
        }
    }
    for (ChartInfo& chart : raw_block)
        for (NoteInfo& timestamp : chart.notes) {
            int next_bpm_index = 0;
            float last_bpm = 120.f;
            float last_bpm_time = 0.f;
            float last_bpm_beat = 0.f;
            while (next_bpm_index < bpms.size() && bpms[next_bpm_index].beat <= timestamp.rowTime) {
                last_bpm_time += (bpms[next_bpm_index].beat - last_bpm_beat) * 240.f / last_bpm;
                last_bpm_beat = bpms[next_bpm_index].beat;
                last_bpm = bpms[next_bpm_index].bpm;
                next_bpm_index += 1;
            }
            timestamp.rowTime = last_bpm_time + (timestamp.rowTime - last_bpm_beat) * 240.f / last_bpm;
        }
    return raw_block;
}

vector<NoteInfo> parse_main_block(stringstream& sm_text) {
    vector<NoteInfo> output;
    int notes;
    int column_value;
    float measure_size = 0.f;
    float measure_number = 0.f;
    vector<int> measure;
    for (std::string line; std::getline(sm_text, line); )
    {
        notes = 0;
        column_value = 1;
        if (line[0] == ',') {
            float inside = 0.f;
            for(unsigned int note_row : measure) {
                if (note_row != 0) {
                    output.push_back(NoteInfo {note_row, measure_number + inside / measure_size});
                }
                inside += 1.f;
            }
            measure_number += 1.f;
            measure.clear();
            measure_size = 0.f;
            continue;
        }
        for(char & it : line) {
            if (it == '\n') {
                break;
            }
            if (isspace(it)) {
                continue;
            }
            if (it == '1' || it == '2') {
                notes += column_value;
            }
            column_value *= 2;
        }
        measure_size += 1.f;
        measure.push_back(notes);
    }
    float inside = 0.f;
    for(unsigned int note_row : measure) {
        if (note_row != 0) {
            output.push_back(NoteInfo{note_row, measure_number + inside / measure_size});
        }
        inside += 1.f;
    }
    measure.clear();

    return output;
}

BPMs parse_bpms_block(stringstream& bpms_block) {
    float next_time;
    float next_bpm;
    BPMs bpm_list = vector<BPM>();
    while (bpms_block.rdbuf()->in_avail()) {
        while (!(bpms_block >> next_time)) {
            if (!bpms_block.rdbuf()->in_avail())
                break;
            bpms_block.get();
        }
        while (bpms_block.get() != '=') {
            if (!bpms_block.rdbuf()->in_avail())
                break;
        }
        while (!(bpms_block >> next_bpm)) {
            if (!bpms_block.rdbuf()->in_avail())
                break;
            bpms_block.get();
        }
        while (bpms_block.get() != ',') {
            if (!bpms_block.rdbuf()->in_avail())
                break;
        }
        bpm_list.push_back(BPM{next_time / 4.f, next_bpm});
    }
    return bpm_list;
}