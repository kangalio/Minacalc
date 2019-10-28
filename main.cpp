#include "minacalc.h"
#include <iostream>

int main() {
    std::vector<NoteInfo> jeff{NoteInfo {3, 0.0}, NoteInfo {3, 40000000.0}, NoteInfo {3, 40000001.0}};
    std::vector<float> x = MinaSDCalc(jeff, 4, 1.f, 0.93f, 1.f);
    std::cout << "Hello, World! this garbage is: " << x[0] << std::endl;
    return 0;
}