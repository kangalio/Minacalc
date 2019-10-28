#include "minacalc.h"
#include <iostream>

int main() {
    std::vector<NoteInfo> jeff;
    jeff.push_back(NoteInfo {3, 0.0});
    jeff.push_back(NoteInfo {3, 40000000.0});
    jeff.push_back(NoteInfo {3, 40000001.0});
    std::cout << "Hello, World! this garbage is: " << MinaSDCalc(jeff, 4, 1.f, 0.93f, 1.f)[0] << std::endl;
    return 0;
}