#include "minacalc.h"
#include "smloader.h"
#include "solocalc.h"
#include <iostream>

using std::cout;
using std::endl;

void printDifficulty(const ChartRating& rating) {
    cout << rating.difficultyName << ":" << endl;
    cout << "Overall: " << rating.rating.overall << endl;
    cout << "Stream: " << rating.rating.stream << endl;
    cout << "JumpStream: " << rating.rating.jumpstream << endl;
    cout << "HandStream: " << rating.rating.handstream << endl;
    cout << "Stamina: " << rating.rating.stamina << endl;
    cout << "Jackspeed: " << rating.rating.jack << endl;
    cout << "Chordjack: " << rating.rating.chordjack << endl;
    cout << "Technical: " << rating.rating.technical;
}

std::vector<ChartRating> difficultyFromFile(const std::string& location) {
    std::ifstream sm_file;
    sm_file.open(location);
    std::vector<ChartRating> rating;
    if (sm_file.is_open()) {
        SMNotes chart = load_from_file(sm_file);
        for (auto& difficulty : chart) {
            rating.push_back(ChartRating { difficulty.difficultyName, MinaSDCalc(difficulty.notes, 1.f, 0.93f)});
        }
    }
    else {
        std::cerr << "failed to open the file" << endl;
    }
    sm_file.close();
    return rating;
}

int main(int argc, char *argv[]) {
    std::vector<ChartRating> rating;
    if (argc > 2) {
        cout << "Solo Difficulty: ";
        std::ifstream sm_file;
        sm_file.open(argv[1]);
        if (sm_file.is_open()) {
            SMNotes chart = load_from_file(sm_file);
            for (auto& difficulty : chart) {
                cout << soloCalc(difficulty.notes, 1.0f, 0.93f) << endl;
            }
        }
    } else if (argc > 1)
        rating = difficultyFromFile(argv[1]);
    else
        rating = difficultyFromFile("../chart.sm");
    //~ for (auto& diff : rating) {
        //~ printDifficulty(diff);
        //~ cout << endl << endl;
    //~ }
    cout << endl;
    cout << "Result: " << rating[rating.size() - 1].overall << endl;
    cout << "Should be: 17.7187" << endl;
    return 0;
}
