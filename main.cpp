#include "minacalc.h"
#include "smloader.h"
#include <iostream>

using std::cout;
using std::endl;

int main() {
    std::ifstream sm_file;
    sm_file.open("../chart.sm");
    if (sm_file.is_open()) {
        std::vector<std::vector<NoteInfo> > chart = load_from_file(sm_file);
        for (auto& difficulty : chart) {
            DifficultyRating rating = MinaSDCalc(difficulty, 1.f, 0.93f);
            cout << "Overall: " << rating.overall << endl;
            cout << "Stream: " << rating.stream << endl;
            cout << "JumpStream: " << rating.jumpstream << endl;
            cout << "HandStream: " << rating.handstream << endl;
            cout << "Stamina: " << rating.stamina << endl;
            cout << "Jackspeed: " << rating.jack << endl;
            cout << "Chordjack: " << rating.chordjack << endl;
            cout << "Technical: " << rating.technical << endl << endl;
        }
    }
    else {
        cout << "failed to open the file" << endl;
    }
    sm_file.close();
    return 0;
}