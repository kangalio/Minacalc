#include "minacalc.h"
#include "smloader.h"
#include <iostream>

using std::cout;
using std::endl;

std::vector<DifficultyRating> difficultyFromFile(const std::string& location) {
    std::ifstream sm_file;
    sm_file.open(location);
    std::vector<DifficultyRating> rating;
    if (sm_file.is_open()) {
        std::vector<std::vector<NoteInfo> > chart = load_from_file(sm_file);
        for (auto& difficulty : chart) {
            rating.push_back(MinaSDCalc(difficulty, 1.f, 0.93f));
            cout << "Overall: " << rating.back().overall << endl;
            cout << "Stream: " << rating.back().stream << endl;
            cout << "JumpStream: " << rating.back().jumpstream << endl;
            cout << "HandStream: " << rating.back().handstream << endl;
            cout << "Stamina: " << rating.back().stamina << endl;
            cout << "Jackspeed: " << rating.back().jack << endl;
            cout << "Chordjack: " << rating.back().chordjack << endl;
            cout << "Technical: " << rating.back().technical << endl << endl;
        }
    }
    else {
        std::cerr << "failed to open the file" << endl;
    }
    sm_file.close();
    return rating;
}

int main(int argc, char *argv[]) {
    if (argc > 1)
        difficultyFromFile(argv[1]);
    else
        difficultyFromFile("../chart.sm");
    return 0;
}