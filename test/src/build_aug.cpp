#include <iostream>

#define VERBOSE

#include <common.hpp>

#include <sdsl/io.hpp>

#include <aug_phoni.hpp>
#include <thr_lce_ds.hpp>

#include <malloc_count.h>

typedef std::pair<std::string, std::vector<uint8_t>> pattern_t;


int main(int argc, char *const argv[]) {
  Args args;
  parseArgs(argc, argv, args);

  verbose("Building the augmented phoni index");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  if (args.thr_lce == "compressed")
  {
    verbose("Using compressed threshold LCEs with byte cap:", (args.bytes) ? std::to_string(args.bytes) : "none");
    ms_pointers<PlainSlp<var_t, Fblc, Fblc>, ms_rle_string_sd, thr_lce_bv<>> ms;
    ms.build(args.filename, args.bytes);

    std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("Augmented PHONI index construction complete");
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    {
    ofstream outfile(args.filename + ".aug", std::ios::binary);
    ms.serialize(outfile);
    }
  }
  else {
    verbose("Using plain threshold LCEs with byte cap:", (args.bytes) ? std::to_string(args.bytes) : "none");
    ms_pointers<> ms;
    ms.build(args.filename, args.bytes);

    std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("Augmented PHONI index construction complete");
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    {
    ofstream outfile(args.filename + ".aug", std::ios::binary);
    ms.serialize(outfile);
    }
  }

  return 0;
}
