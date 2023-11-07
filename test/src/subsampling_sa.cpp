#include "sdsl/int_vector.hpp"
#include <cstdint>
#include <iostream>
#include <ostream>
#include <vector>

#define VERBOSE

#include <common.hpp>
#include <subsampling_sa.hpp>

uint64_t performSubSampling(const size_t &maxLF, const sdsl::int_vector<> &samples,
                        ri::sparse_sd_vector &subsampled_runIdx_bv_sparse,
                        sdsl::int_vector<> &subsampled_values);

int main(int argc, char *const argv[]) {
  Args args;
  parseArgs(argc, argv, args);

  // Subsample the samples of the run starts
  {
    verbose("Subsampling the samples of the run starts");
    std::chrono::high_resolution_clock::time_point t_insert_start =
        std::chrono::high_resolution_clock::now();

    sdsl::int_vector<> samples_start;
    read_samples(args.filename + ".ssa", samples_start);
    ri::sparse_sd_vector subsampled_starts_runIdx_bv_sparse;
    sdsl::int_vector<> subsampled_starts_values;
    uint64_t ssa_r = performSubSampling(args.maxLF, samples_start,
                       subsampled_starts_runIdx_bv_sparse,
                       subsampled_starts_values);
    std::cout << "ssa_samples: " << std::to_string(ssa_r) << std::endl;

    string outputFileStartValues = args.filename + ".s." + std::to_string(args.maxLF) + ".ssa";
    std::ofstream outStartValues(outputFileStartValues, std::ios::binary);
    subsampled_starts_values.serialize(outStartValues);

    string outputFileStartBv = args.filename + ".s." + std::to_string(args.maxLF) + ".ssa.bv";
    std::ofstream outStartBv(outputFileStartBv, std::ios::binary);
    subsampled_starts_runIdx_bv_sparse.serialize(outStartBv);

    std::chrono::high_resolution_clock::time_point t_insert_end =
        std::chrono::high_resolution_clock::now();

    verbose("samples_start subsampling complete");
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(
                                      t_insert_end - t_insert_start)
                                      .count());
  }

  // Subsample the samples of the run ends

  {
    verbose("Subsampling the samples of the run ends");
    std::chrono::high_resolution_clock::time_point t_insert_start =
        std::chrono::high_resolution_clock::now();

    sdsl::int_vector<> samples_end;
    read_samples(args.filename + ".esa", samples_end);
    ri::sparse_sd_vector subsampled_ends_runIdx_bv_sparse;
    sdsl::int_vector<> subsampled_ends_values;
    uint64_t esa_r = performSubSampling(args.maxLF, samples_end,
                       subsampled_ends_runIdx_bv_sparse,
                       subsampled_ends_values);
    std::cout << "esa_samples: " << std::to_string(esa_r) << std::endl;

    string outputFileEndValues = args.filename + ".s." + std::to_string(args.maxLF) + ".esa";
    std::ofstream outEndValues(outputFileEndValues, std::ios::binary);
    subsampled_ends_values.serialize(outEndValues);

    string outputFileEndBv = args.filename + ".s." + std::to_string(args.maxLF) + ".esa.bv";
    std::ofstream outEndBv(outputFileEndBv, std::ios::binary);
    subsampled_ends_runIdx_bv_sparse.serialize(outEndBv);

    std::chrono::high_resolution_clock::time_point t_insert_end =
        std::chrono::high_resolution_clock::now();

    verbose("samples_last subsampling complete");
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(
                                      t_insert_end - t_insert_start)
                                      .count());
  }
}

uint64_t performSubSampling(const size_t &maxLF, const sdsl::int_vector<> &samples,
                        ri::sparse_sd_vector &subsampled_runIdx_bv_sparse,
                        sdsl::int_vector<> &subsampled_values) {

  uint64_t r = samples.size();

  std::vector<uint64_t> samples_runIdx_sorted;
  samples_runIdx_sorted.resize(r);

  for (uint64_t i = 0; i < r; i++) {
    samples_runIdx_sorted[i] = i;
  }

  sort(samples_runIdx_sorted.begin(), samples_runIdx_sorted.end(),
       [&samples](const uint64_t &a, const uint64_t &b) -> bool {
         return samples[a] < samples[b];
       });

  std::vector<uint64_t> subsampled_runIdxs; // Indices of subsampled BWT runs

  uint64_t r_prime = 0;

  // Compute sampled BWT run tails
  uint64_t last_subsampled_runIdx = samples_runIdx_sorted[0];
  uint64_t last_subsampled_value = samples[last_subsampled_runIdx];
  subsampled_runIdxs.emplace_back(last_subsampled_runIdx);

  for (uint64_t i = 1; i < r; i++) {
    uint64_t current_runIdx = samples_runIdx_sorted[i];
    uint64_t current_value = samples[current_runIdx];

    if (current_value - last_subsampled_value >= maxLF) {
      subsampled_runIdxs.emplace_back(current_runIdx);
      last_subsampled_runIdx = current_runIdx;
      last_subsampled_value = current_value;
    }
  }

  r_prime = subsampled_runIdxs.size();
  verbose("Maximum number of LF steps between two samples: " + std::to_string(maxLF));
  verbose("r = " + std::to_string(r));
  verbose("r_prime = " + std::to_string(r_prime));

  // Sort indexes of subsampled runs
  sort(subsampled_runIdxs.begin(), subsampled_runIdxs.end());

  auto log_n = sdsl::bits::hi(last_subsampled_value) + 1;

  // Text position of subsampled BWT runs in order of the BWT Array
  subsampled_values = sdsl::int_vector<>(r_prime, 0, log_n);
  // Mark which runs are sampled (bitvector)
  sdsl::bit_vector subsampled_runIdx_bv = sdsl::bit_vector(r, 0);

  // Compute text position of the subsampled BWT runs and mark the sampled
  // indices in a bitvector
  for (uint64_t i = 0; i < r_prime; i++) {
    uint64_t runIdx = subsampled_runIdxs[i];
    uint64_t run_value = samples[runIdx];
    subsampled_values[i] = run_value;
    subsampled_runIdx_bv[runIdx] = 1;
  }

  // Mark which runs are sampled (sparse bitvector)
  subsampled_runIdx_bv_sparse = ri::sparse_sd_vector(subsampled_runIdx_bv);
  return r_prime;
}