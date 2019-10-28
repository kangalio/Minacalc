#include "minacalc.h"
#include "smloader.h"
#include <iostream>

int main() {
    ifstream sm_file;
    sm_file.open("../chart.sm");
    if (sm_file.is_open()) {
        std::vector<NoteInfo> jeff = load_from_file(sm_file);
        std::vector<float> x = MinaSDCalc(jeff, 4, 1.f, 0.93f, 1.f);
        std::cout << "Overall: " << x[0] << std::endl;
        std::cout << "Stream: " << x[1] << std::endl;
        std::cout << "JumpStream: " << x[2] << std::endl;
        std::cout << "HandStream: " << x[3] << std::endl;
        std::cout << "Stamina: " << x[4] << std::endl;
        std::cout << "Jackspeed: " << x[5] << std::endl;
        std::cout << "Chordjack: " << x[6] << std::endl;
        std::cout << "Technical: " << x[7] << std::endl;
    }
    else {
        cout << "failed to open the file" << std::endl;
    }
    sm_file.close();
    return 0;
}