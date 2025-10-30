#include <iostream>
#include <random>

struct Ising {
  int L;
  std::vector<int> spins;
  std::mt19937 rng;
  std::uniform_real_distribution<double> U01;

  Ising(int L, unsigned seed = 67)
      : L(L), spins(L * L), rng(seed), U01(0.0, 1.0) {
    for (int i = 0; i < L * L; i++) {
      spins[i] = (U01(rng) < 0.5) ? 1 : -1;
    }
  }

  // Gets the index associated with these coordinates, wrapping if needed
  inline int get_idx(int x, int y) { return (x + L) % L + ((y + L) % L) * L; }
};

int main(int argc, char *argv[]) {
  Ising ising(5);

  for (int i = 0; i < ising.L; i++) {
    for (int j = 0; j < ising.L; j++) {
      std::cout << ((ising.spins[i * ising.L + j] == 1) ? "+" : "-") << " ";
    }
    std::cout << "\n";
  }

  return 0;
}
