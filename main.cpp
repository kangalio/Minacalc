#include "minacalc.h"
#include "smloader.h"
#include <iostream>

int main() {
    ifstream sm_file;
    sm_file.open("../chart.sm");
    if (sm_file.is_open()) {
        std::vector<NoteInfo> jeff = load_from_file(sm_file);
        for (NoteInfo note : jeff) {
            cout << note.notes << " " << note.rowTime << endl;
        }
        DifficultyRating x = MinaSDCalc(jeff, 1.f, 0.93f);
        std::cout << "Overall: " << x.overall << std::endl;
        std::cout << "Stream: " << x.stream << std::endl;
        std::cout << "JumpStream: " << x.jumpstream << std::endl;
        std::cout << "HandStream: " << x.handstream << std::endl;
        std::cout << "Stamina: " << x.stamina << std::endl;
        std::cout << "Jackspeed: " << x.jack << std::endl;
        std::cout << "Chordjack: " << x.chordjack << std::endl;
        std::cout << "Technical: " << x.technical << std::endl;
    }
    else {
        cout << "failed to open the file" << std::endl;
    }
    sm_file.close();
    return 0;
}