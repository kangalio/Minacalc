#include "minacalc.h"
#include "smloader.h"
#include <iostream>

using std::cout;
using std::endl;

int main() {
    ifstream sm_file;
    sm_file.open("../chart.sm");
    if (sm_file.is_open()) {
        std::vector<NoteInfo> jeff = load_from_file(sm_file);
        for (NoteInfo note : jeff) {
            cout << note.notes << " " << note.rowTime << endl;
        }
        DifficultyRating x = MinaSDCalc(jeff, 1.f, 0.93f);
        cout << "Overall: " << x.overall << endl;
        cout << "Stream: " << x.stream << endl;
        cout << "JumpStream: " << x.jumpstream << endl;
        cout << "HandStream: " << x.handstream << endl;
        cout << "Stamina: " << x.stamina << endl;
        cout << "Jackspeed: " << x.jack << endl;
        cout << "Chordjack: " << x.chordjack << endl;
        cout << "Technical: " << x.technical << endl;
    }
    else {
        cout << "failed to open the file" << endl;
    }
    sm_file.close();
    return 0;
}