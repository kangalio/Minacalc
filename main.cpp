#include "minacalc.h"
#include "smloader.h"
#include "solocalc.h"
#include <iostream>

using std::cout;
using std::endl;

void printDifficulty(const DifficultyRating& rating) {
    cout << "Overall: " << rating.overall << endl;
    cout << "Stream: " << rating.stream << endl;
    cout << "JumpStream: " << rating.jumpstream << endl;
    cout << "HandStream: " << rating.handstream << endl;
    cout << "Stamina: " << rating.stamina << endl;
    cout << "Jackspeed: " << rating.jack << endl;
    cout << "Chordjack: " << rating.chordjack << endl;
    cout << "Technical: " << rating.technical;
}

std::vector<DifficultyRating> difficultyFromFile(const std::string& location) {
    std::ifstream sm_file;
    sm_file.open(location);
    std::vector<DifficultyRating> rating;
    if (sm_file.is_open()) {
        std::vector<std::vector<NoteInfo> > chart = load_from_file(sm_file);
        for (auto& difficulty : chart) {
            rating.push_back(MinaSDCalc(difficulty, 1.f, 0.93f));
        }
    }
    else {
        std::cerr << "failed to open the file" << endl;
    }
    sm_file.close();
    return rating;
}

int main(int argc, char *argv[]) {
    std::vector<DifficultyRating> rating;
    if (argc > 2) {
        cout << "Solo Difficulty: ";
        std::ifstream sm_file;
        sm_file.open(argv[1]);
        if (sm_file.is_open()) {
            std::vector<std::vector<NoteInfo> > chart = load_from_file(sm_file);
            for (auto& difficulty : chart) {
                cout << soloCalc(difficulty, 1.0f, 0.93f) << endl;
            }
        }
    } else if (argc > 1)
        rating = difficultyFromFile(argv[1]);
    else
        rating = difficultyFromFile("../chart.sm");
    for (auto& diff : rating) {
        printDifficulty(diff);
        cout << endl << endl;
    }
    return 0;
}