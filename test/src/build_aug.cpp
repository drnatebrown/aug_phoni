#include <iostream>

#define VERBOSE

#include <common.hpp>

#include <sdsl/io.hpp>

#include <lazy_lce.hpp>

#include <malloc_count.h>

typedef std::pair<std::string, std::vector<uint8_t>> pattern_t;


int main(int argc, char *const argv[]) {
  Args args;
  parseArgs(argc, argv, args);

  verbose("Building the augmented phoni index");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  verbose("Using plain threshold LCEs with byte cap:", (args.bytes) ? std::to_string(args.bytes) : "none");
  ms_pointers<> ms;
  ms.build(args.filename, args.maxLF, args.bytes);

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Augmented PHONI index construction complete");
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  {
  ofstream outfile(args.filename + "." + to_string(args.maxLF) + ".aug.mlq", std::ios::binary);
  ms.serialize(outfile);
  }

  return 0;
}
