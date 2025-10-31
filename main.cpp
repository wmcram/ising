#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <thread>

struct Ising {
  int L;
  double T;
  double J;
  std::vector<int> spins;
  std::mt19937 rng;
  std::uniform_real_distribution<double> U01;

  Ising(int L, double T, double J,
        unsigned seed = std::chrono::high_resolution_clock::now()
                            .time_since_epoch()
                            .count())
      : L(L), T(T), J(J), spins(L * L), rng(seed), U01(0.0, 1.0) {
    for (int i = 0; i < L * L; i++) {
      spins[i] = (U01(rng) < 0.5) ? 1 : -1;
    }
  }

  // Gets the index associated with these coordinates, wrapping if needed
  inline int get_idx(int x, int y) const {
    return (x + L) % L + ((y + L) % L) * L;
  }

  // Gets the energy for a single site (x,y)
  int neighbor_sum(int x, int y) const {
    return spins[get_idx(x + 1, y)] + spins[get_idx(x - 1, y)] +
           spins[get_idx(x, y + 1)] + spins[get_idx(x, y - 1)];
  }

  // Gets the total energy of the system
  double total_energy() const {
    double E = 0.0;
    for (int x = 0; x < L; x++) {
      for (int y = 0; y < L; y++) {
        int si = spins[get_idx(x, y)];
        int su = spins[get_idx(x, y + 1)];
        int sr = spins[get_idx(x + 1, y)];

        E -= J * si * sr;
        E -= J * si * su;
      }
    }
    return E;
  }

  // Performs a metropolis sweep, updating the spins
  void sweep() {
    double beta = 1.0 / T;
    std::uniform_int_distribution<int> pick(0, L - 1);
    for (int n = 0; n < L * L; n++) {
      int x = pick(rng);
      int y = pick(rng);
      int idx = get_idx(x, y);
      int si = spins[idx];
      int nbr = neighbor_sum(x, y);
      double dE = 2 * si * J * nbr;
      if (dE <= 0.0) {
        spins[idx] = -si;
      } else {
        double r = U01(rng);
        if (r < exp(-beta * dE))
          spins[idx] = -si;
      }
    }
  }

  // Prints the state of the spins to the terminal.
  void print_state() {
    std::cout << "\x1B[2J\x1B[H";
    for (int i = 0; i < L; i++) {
      for (int j = 0; j < L; j++) {
        std::cout << ((spins[get_idx(i, j)] == 1) ? "+" : "-") << " ";
      }
      std::cout << "\n";
    }
    std::cout << "T: " << T << " E: " << total_energy() << std::endl;
  }

  // Iterates the sweep() method for a number of rounds, printing the state and
  // then sleeping in between each round for the specified delay.
  void sweep_iter(int rounds, int delay = 0, double delta_temp = 0.0) {
    for (int i = 0; i < rounds; i++) {
      sweep();
      print_state();
      if (delay > 0) {
        std::this_thread::sleep_for(std::chrono::milliseconds(delay));
      }
      T += delta_temp;
    }
  }
};

int main(int argc, char *argv[]) {
  Ising ising(27, 2.5, 1.0);
  ising.sweep_iter(1000, 250, -0.01);
  return 0;
}
