// Link against the release LIBOBJ files (excluding main.o and split.o).
#include "vina.h"
#include <chrono>
#include <cstdlib>
#include <iostream>

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: benchmark_ligand_loading ligand.pdbqt iterations vina|vinardo|ad4 single|vector\n";
        return 2;
    }
    const int iterations = std::atoi(argv[2]);
    const std::string mode(argv[4]);
    if (iterations < 1 || (mode != "single" && mode != "vector")) return 2;
    Vina vina(argv[3], 1, 42, 0);
    const std::string ligand = get_file_contents(argv[1]);
    const std::vector<std::string> ligands(1, ligand);
    // Warm-up and replacement of an already populated table are both exercised.
    if (mode == "single") vina.set_ligand_from_string(ligand);
    else vina.set_ligand_from_string(ligands);
    const auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < iterations; ++i) {
        if (mode == "single") vina.set_ligand_from_string(ligand);
        else vina.set_ligand_from_string(ligands);
    }
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    std::cout << seconds / iterations << '\n';
}
