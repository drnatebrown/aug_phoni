/* ms_pointers - Computes the matching statistics pointers from BWT and Thresholds
    Copyright (C) 2020 Massimiliano Rossi

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*!
   \file ms_pointers.hpp
   \brief ms_pointers.hpp Computes the matching statistics pointers from BWT and Thresholds.
   \author Massimiliano Rossi
   \date 09/07/2020
*/

#ifndef _AUG_MS_POINTERS_HH
#define _AUG_MS_POINTERS_HH


/** FLAGS **/
#define MEASURE_TIME 1  //measure the time for LCE and backward search?
//#define NAIVE_LCE_SCHEDULE 1 //stupidly execute two LCEs without heurstics
#include "Common.hpp"
#include <cstddef>
#include <sdsl/util.hpp>
#ifndef NAIVE_LCE_SCHEDULE //apply a heuristic
#define SORT_BY_DISTANCE_HEURISTIC 1 // apply a heuristic to compute the LCE with the closer BWT position first
#endif

#ifndef DCHECK_HPP
#define DCHECK_HPP
#include <string>
#include <sstream>
#include <stdexcept>

#ifndef DCHECK
#ifdef NDEBUG
#define ON_DEBUG(x)
#define DCHECK_(x,y,z)
#define DCHECK(x)
#define DCHECK_EQ(x, y)
#define DCHECK_NE(x, y)
#define DCHECK_LE(x, y)
#define DCHECK_LT(x, y)
#define DCHECK_GE(x, y)
#define DCHECK_GT(x, y)
#else//NDEBUG
#define ON_DEBUG(x) x
#define DCHECK_(x,y,z) \
  if (!(x)) throw std::runtime_error(std::string(" in file ") + __FILE__ + ':' + std::to_string(__LINE__) + (" the check failed: " #x) + ", we got " + std::to_string(y) + " vs " + std::to_string(z))
#define DCHECK(x) \
  if (!(x)) throw std::runtime_error(std::string(" in file ") + __FILE__ + ':' + std::to_string(__LINE__) + (" the check failed: " #x))
#define DCHECK_EQ(x, y) DCHECK_((x) == (y), x,y)
#define DCHECK_NE(x, y) DCHECK_((x) != (y), x,y)
#define DCHECK_LE(x, y) DCHECK_((x) <= (y), x,y)
#define DCHECK_LT(x, y) DCHECK_((x) < (y) ,x,y)
#define DCHECK_GE(x, y) DCHECK_((x) >= (y),x,y )
#define DCHECK_GT(x, y) DCHECK_((x) > (y) ,x,y)
#endif //NDEBUG
#endif //DCHECK
#endif /* DCHECK_HPP */


#include <common.hpp>

#include <malloc_count.h>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>

#include <r_index.hpp>

#include <ms_rle_string.hpp>
#include <thresholds_ds.hpp>
#include <thr_lce_ds.hpp>

#include "PlainSlp.hpp"
#include "PoSlp.hpp"
#include "ShapedSlp_Status.hpp"
#include "ShapedSlp.hpp"
#include "ShapedSlpV2.hpp"
#include "SelfShapedSlp.hpp"
#include "SelfShapedSlpV2.hpp"
#include "DirectAccessibleGammaCode.hpp"
#include "IncBitLenCode.hpp"
#include "FixedBitLenCode.hpp"
#include "SelectType.hpp"
#include "VlcVec.hpp"

#ifdef MEASURE_TIME
struct Stopwatch {
    std::chrono::high_resolution_clock::time_point m_val = std::chrono::high_resolution_clock::now();

    void reset()
    {
        m_val = std::chrono::high_resolution_clock::now();
    }
    double seconds() const
    {
        return std::chrono::duration<double, std::ratio<1>>(std::chrono::high_resolution_clock::now() - m_val).count();
    }

};
#endif//MEASURE_TIME

using var_t = uint32_t;
using Fblc = FixedBitLenCode<>;
using SelSd = SelectSdvec<>;
using SelMcl = SelectMcl<>;
using DagcSd = DirectAccessibleGammaCode<SelSd>;
using DagcMcl = DirectAccessibleGammaCode<SelMcl>;
using Vlc64 = VlcVec<sdsl::coder::elias_delta, 64>;
using Vlc128 = VlcVec<sdsl::coder::elias_delta, 128>;


template <
    class SlpT = PlainSlp<var_t, Fblc, Fblc>,
    class rle_string_t = ms_rle_string_sd,
    class thr_lce_t = thr_lce_plain<rle_string_t>,
    class sparse_bv_type = ri::sparse_sd_vector,
    class thresholds_t = thr_bv<rle_string_t>
    >
class ms_pointers : ri::r_index<sparse_bv_type, rle_string_t>
{
public:
    thresholds_t thresholds;
    thr_lce_t thr_lce;

    SlpT slp;

    int_vector<> samples_start;
    std::queue<ulint> LF_next;
    sparse_bv_type subsampled_start_samples_bv;
    sparse_bv_type subsampled_last_samples_bv;

    typedef size_t size_type;

    ms_pointers()
        : ri::r_index<sparse_bv_type, rle_string_t>()
    {}

    void build(const std::string& filename, int max_LF = 0, int bytes = 0)
    {
        verbose("Building the r-index from BWT");

        std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

        std::string bwt_fname = filename + ".bwt";

        verbose("RLE encoding BWT and computing SA samples");

        if (true) {
            std::string bwt_heads_fname = bwt_fname + ".heads";
            std::ifstream ifs_heads(bwt_heads_fname);
            std::string bwt_len_fname = bwt_fname + ".len";
            std::ifstream ifs_len(bwt_len_fname);
            this->bwt = rle_string_t(ifs_heads, ifs_len);

            ifs_heads.seekg(0);
            ifs_len.seekg(0);
            this->build_F_(ifs_heads, ifs_len);
        } else {
            std::ifstream ifs(bwt_fname);
            this->bwt = rle_string_t(ifs);

            ifs.seekg(0);
            this->build_F(ifs);
        }

        this->r = this->bwt.number_of_runs();
        ri::ulint n = this->bwt.size();
        int log_n = bitsize(uint64_t(this->bwt.size()));

        verbose("Number of BWT equal-letter runs: r = ", this->r);
        verbose("Rate n/r = ", double(this->bwt.size()) / this->r);
        verbose("log2(r) = ", log2(double(this->r)));
        verbose("log2(n/r) = ", log2(double(this->bwt.size()) / this->r));
        verbose("log2(n) = ", log_n);

        string subsamples_start_filename = filename + ".s." + std::to_string(max_LF) + ".ssa";
        ifstream subsamples_start_file(subsamples_start_filename, ios::binary);
        samples_start.load(subsamples_start_file);

        string subsamples_start_bv_filename = filename + ".s." + std::to_string(max_LF) + ".ssa.bv";
        ifstream subsamples_start_bv_file(subsamples_start_bv_filename, ios::binary);
        subsampled_start_samples_bv.load(subsamples_start_bv_file);

        string subsamples_last_filename = filename + ".s." + std::to_string(max_LF) + ".esa";
        ifstream subsamples_last_file(subsamples_last_filename, ios::binary);
        this->samples_last.load(subsamples_last_file);

        string subsamples_last_bv_filename = filename + ".s." + std::to_string(max_LF) + ".esa.bv";
        ifstream subsamples_last_bv_file(subsamples_last_bv_filename, ios::binary);
        subsampled_last_samples_bv.load(subsamples_last_bv_file);

        std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();
        verbose("R-index construction complete");
        verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
        verbose(3);

        //load_grammar(filename);

        verbose("text length: ", slp.getLen());
        verbose("bwt length: ", this->bwt.size());

        thresholds = thresholds_t(filename, &this->bwt);
        thr_lce = thr_lce_t(filename, &this->bwt, bytes);

        verbose("finished augmented thresholds construction");
    }


    void load_grammar(const std::string& filename)
    {
        {
            verbose("Load Grammar");
            std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();
            if (std::is_same<SlpT, PlainSlp<var_t, Fblc, Fblc>>::value) {
                ifstream fs(filename + ".plain.slp");
                slp.load(fs);
                fs.close();
            } else {
                ifstream fs(filename + ".slp");
                slp.load(fs);
                fs.close();
            }

            std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();
            verbose("Memory peak: ", malloc_count_peak());
            verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
        }
        DCHECK_EQ(slp.getLen()+1, this->bwt.size());
    }



    void read_samples(std::string filename, ulint r, ulint n, int_vector<> &samples)
    {
        int log_n = bitsize(uint64_t(n));
        struct stat filestat;
        FILE *fd;

        if ((fd = fopen(filename.c_str(), "r")) == nullptr)
            error("open() file " + filename + " failed");

        int fn = fileno(fd);
        if (fstat(fn, &filestat) < 0)
            error("stat() file " + filename + " failed");

        if (filestat.st_size % SSABYTES != 0)
            error("invilid file " + filename);

        size_t length = filestat.st_size / (2*SSABYTES);
        //Check that the length of the file is 2*r elements of 5 bytes
        assert(length == r);

        // Create the vector
        samples = int_vector<>(r, 0, log_n);

        // Read the vector
        uint64_t left = 0;
        uint64_t right = 0;
        size_t i = 0;
        while (fread((char *)&left, SSABYTES, 1, fd) && fread((char *)&right, SSABYTES, 1,fd)) {
            ulint val = (right ? right - 1 : n - 1);
            assert(bitsize(uint64_t(val)) <= log_n);
            samples[i++] = val;
        }

        fclose(fd);
    }

    void write_int(ostream& os, const size_t& i)
    {
        os.write(reinterpret_cast<const char*>(&i), sizeof(size_t));
    }

    // Computes the matching statistics pointers for the given pattern
    size_t query(const std::string& patternfile, const std::string& len_filename, const std::string& ref_filename) {
        const char* p;
        size_t m;
        map_file(patternfile.c_str(), p, m);

        auto pattern_at = [&] (const size_t pos) {
            return p[pos];
        };

        ofstream len_file(len_filename, std::ios::binary);
        ofstream ref_file(ref_filename, std::ios::binary);

        const size_t n = slp.getLen();
        verbose("pattern length: ", m);

        size_t last_len;
        size_t last_ref;

        size_t lce_freq = 3, lce_cnt = 0; // freq of LCE skipping and a counter for it
        vector<size_t> stored_sample_pos(lce_freq), stored_ptr(lce_freq);
        vector<int> stored_it(lce_freq);
        vector<bool> direction(lce_freq);
        bool lce_is_paused = false;

        auto write_len = [&] (const size_t i, bool lce_is_paused=false) {
            if (lce_is_paused) return;
            last_len = i;
            write_int(len_file, last_len);
        };
        auto write_ref = [&] (const size_t i) {
            last_ref = i;
            write_int(ref_file, last_ref);
        };

        write_len(1, lce_is_paused);

        //! Start with the last character Adrian: do we assume the last char occurs in T??
        auto pos = this->bwt.select(1, pattern_at(m-1));
        {
            const ri::ulint run_of_j = this->bwt.run_of_position(pos);
            write_ref(samples_start[run_of_j]);
        }
        pos = LF(pos, pattern_at(m-1));

        #ifdef MEASURE_TIME
        double time_lce = 0;
        double time_backwardstep = 0;
        size_t count_lce_total = 0;
        size_t count_lce_skips = 0;
        #endif

        for (size_t i = 1; i < m; ++i) {
            //verbose("i= ", i);

            const auto c = pattern_at(m - i - 1);
            const size_t number_of_runs_of_c = this->bwt.number_of_letter(c);
            if (number_of_runs_of_c == 0) {
                write_len(0, lce_is_paused);
                write_ref(1); // why 1???
            } else if (pos < this->bwt.size() && this->bwt[pos] == c) {
                write_len(last_len+1, lce_is_paused);
                write_ref(last_ref-1);
            } else { // we jump
                const ri::ulint rank = this->bwt.rank(pos, c);

                ri::ulint run0, run1;
                size_t sa0, sa1;

                if (rank < number_of_runs_of_c) {
                    sa1 = this->bwt.select(rank, c);
                }
                if (rank >= number_of_runs_of_c) {
                    sa0 = this->bwt.select(rank-1, c);
                }

                struct Triplet {
                    size_t sa, ref, len;
                };

                auto delay_succeeding_lce = [&] () -> Triplet {
                    const size_t textposStart = this->samples_start[run1];

                    stored_it[lce_cnt] = i;
                    stored_sample_pos[lce_cnt] = textposStart;
                    stored_ptr[lce_cnt] = last_ref;
                    direction[lce_cnt] = 1;
                    lce_cnt++;
                    return {sa1, textposStart, 0};

                };

                auto delay_preceding_lce = [&] () -> Triplet {
                    const ri::ulint run0 = this->bwt.run_of_position(sa0);
                    const size_t textposLast = this->samples_last[run0];

                    stored_it[lce_cnt] = i;
                    stored_sample_pos[lce_cnt] = textposLast;
                    stored_ptr[lce_cnt] = last_ref;
                    direction[lce_cnt] = 0;
                    lce_cnt++;
                    return {sa0, textposLast, 0};
                };

                auto compute_succeeding_lce = [&] (const size_t textposStart, size_t last_ref, size_t last_len) -> size_t {
                    const size_t lenStart = textposStart+1 >= n ? 0 : lceToRBounded(slp, textposStart+1, last_ref, last_len);
                    const size_t len1 = std::min(last_len, lenStart);
                    return len1;
                };

                auto compute_preceding_lce = [&] (const size_t textposLast, size_t last_ref, size_t last_len) -> size_t {
                    const size_t lenLast = textposLast+1 >= n ? 0 : lceToRBounded(slp, textposLast+1, last_ref, last_len);
                    const size_t len0 = std::min(last_len, lenLast);
                    return len0;
                };

                const Triplet jump = [&] () -> Triplet {
                    if(rank == 0) {
                        // Check only succeeding -> we ignore thresholds in this case
                        return delay_succeeding_lce();
                    } else if(rank >= number_of_runs_of_c) {
                        // Check only preceding -> we ignore thresholds in this case
                        return delay_preceding_lce();
                    }
                    // Check thresholds and boundary LCEs first
                    run1 = this->bwt.run_of_position(sa1);
                    const size_t thr = thresholds[run1];
                    if (pos < thr) {
                        // return delay_preceding_lce();
                        if (!lce_is_paused && thr_lce.skip_preceding_lce(run1, last_len)) {
                            const ri::ulint run0 = this->bwt.run_of_position(sa0);
                            const size_t ref0 = this->samples_last[run0];
                            const size_t len0 = last_len;
                            return {sa0, ref0, len0};
                        } else {
                            return delay_preceding_lce();
                        }
                    } else {
                        // return delay_succeeding_lce();
                        if (!lce_is_paused && thr_lce.skip_succeeding_lce(run1, last_len)) {
                            const size_t ref1 = this->samples_start[run1];;
                            const size_t len1 = last_len;
                            return {sa1, ref1, len1};
                        } else {
                            return delay_succeeding_lce();
                        }
                    }
                }();

                if (lce_cnt == lce_freq) {
                    lce_is_paused = false;

                    size_t lce = 0;
                    if (direction[lce_freq-1] == 0)     lce = compute_preceding_lce(stored_sample_pos[lce_freq-1],last_ref, last_len);
                    else                                lce = compute_succeeding_lce(stored_sample_pos[lce_freq-1], last_ref, last_len);

                    int skipped_steps = (stored_it[lce_freq-1] - stored_it[0]);
                    // if we still process the same MEM, we know the lens right ahead
                    if ((last_len + skipped_steps) == lce) {
                        for (int j = 0; j < skipped_steps; j++) write_len(last_len+1, lce_is_paused); 
                    } else {
                        // do the delayed LCEs
                        for (int j = 0; j < lce_freq; j++) {
                            // TODO: try to use the stored LCE's
                            if (direction[j] == 0)  last_len = compute_preceding_lce(stored_sample_pos[j],last_ref, last_len);
                            else                    last_len = compute_succeeding_lce(stored_sample_pos[j], last_ref, last_len);

                            if (j < (lce_freq - 1)) {
                                for (int k = 0; k < (stored_it[j + 1] - stored_it[j]); k++) {
                                    write_len(last_len+1, lce_is_paused);
                                }
                            }
                        }
                        // the last LCE wasn't written yet
                        write_len(last_len, lce_is_paused);
                        lce_cnt = 0;
                    }
                } else {
                    write_len(1 + jump.len, lce_is_paused);
                }
                write_ref(jump.ref);
                pos = jump.sa;
            }
            pos = LF(pos, c); //! Perform one backward step
        }
#ifdef MEASURE_TIME
        cout << "Time backwardsearch: " << time_backwardstep << std::endl;
        cout << "Time lce: " << time_lce << std::endl;
        cout << "Count lce: " << count_lce_total << std::endl;
        cout << "Count lce skips: " << count_lce_skips << std::endl;
#endif
        return m;
    }
    // Computes the matching statistics pointers for the given pattern
    std::pair<std::vector<size_t>, std::vector<size_t>> query(const std::vector<uint8_t> &pattern)
    {
        size_t m = pattern.size();

        return _query(pattern.data(), m);
    }

    std::pair<std::vector<size_t>, std::vector<size_t>> query(const char *pattern, const size_t m, const size_t lcp_delay)
    {
        return _query(pattern, m, lcp_delay);
    }

    // Computes the matching statistics pointers for the given pattern
    template <typename string_t>
    std::pair<std::vector<size_t>, std::vector<size_t>> _query(const string_t &pattern, const size_t m, const size_t min_ms_len)
    {
        #ifdef MEASURE_TIME
        double time_lce = 0;
        double time_mlq = 0;
		double time_backwardstep = 0;
		size_t count_lce_total = 0;
		size_t count_lce_skips = 0;
        size_t count_mlq_total = 0;
        #endif

        auto pattern_at = [&] (const size_t pos) {
            return pattern[pos];
        };

        std::vector<size_t> lens;
        std::vector<size_t> refs;

        lens.reserve(m);
        refs.reserve(m);

        const size_t n = slp.getLen();
        verbose("Pattern length: ", m);
        slp.precompute_pattern(pattern); // KR hashes of the pattern prefixes

        size_t last_len;
        size_t last_ref;

        size_t skip_queries = 0;
        bool do_mlq = false;

        auto write_len = [&] (const size_t l) {
            last_len = l;
            lens.push_back(last_len);
        };
        auto write_ref = [&] (const size_t p) {
            last_ref = p;
            refs.push_back(last_ref);
        };

        skip_queries = (min_ms_len == 0) ? 0 : min_ms_len - 1;
        if (!skip_queries) write_len(1);

        //! Start with the last character
        auto pos = this->bwt.select(1, pattern_at(m-1));
        {
            if (!skip_queries)
            {
                const ri::ulint run_of_j = this->bwt.run_of_position(pos);
                write_ref(subsamples_start(run_of_j, pos));
            }
        }
        pos = LF(pos, pattern_at(m-1));
        --skip_queries;

        struct Triplet {
            size_t sa, ref, len;
        };

        // auto store_lce_info = [&] (const int i, const size_t textpos, const size_t run, lce_skip skip) {
        //     lce_is_paused = true;
        //     stored_it[lce_cnt] = i;
        //     stored_sample_pos[lce_cnt] = textpos;
        //     stored_ptr[lce_cnt] = last_ref;
        //     stored_run[lce_cnt] = run;
        //     direction[lce_cnt] = skip;
        //     lce_cnt++;
        // };

        // auto delay_preceding_lce = [&] (const size_t run, const size_t rank, char c, int i, lce_skip skip) -> Triplet {
        //     const size_t sa0 = this->bwt.select(rank-1, c);
        //     const ri::ulint run0 = this->bwt.run_of_position(sa0);
        //     const size_t textpos = this->samples_last[run0];
        //     //verbose("i = ", i, "last_ref = ", last_ref, " lce_is_paused =", lce_is_paused, "UP= ", textposLast);
        //     store_lce_info(i, textpos, run, skip);
        //     return {sa0, textpos, 0};
        // };

        // auto delay_succeeding_lce = [&] (const size_t sa, const size_t run, int i, lce_skip skip) -> Triplet {
        //     const size_t textpos = this->samples_start[run];
        //     //verbose("i = ", i, "last_ref = ", last_ref, " lce_is_paused =", lce_is_paused, "DOWN= ", textposStart, " run = ", this->bwt.run_of_position(sa1));
        //     store_lce_info(i, textpos, run, skip);
        //     return {sa, textpos, 0};
        // };

        auto compute_mlq = [&] (const size_t pos_sample, const size_t pos_pattern) -> size_t {
            //verbose("computing MLQ for:  ", pos_sample, " ", pos_pattern);
            if (pos_sample + 1 >= n)
                return 0;
            else
            {
                #ifdef MEASURE_TIME
                Stopwatch s;
                #endif
                const size_t lenLast = match_length_query(slp, pos_sample + 1, pos_pattern);
                #ifdef MEASURE_TIME
                time_mlq += s.seconds();
                ++count_mlq_total;
                #endif
                return lenLast;
            }
        };

        auto compute_lce = [&] (const size_t pos_sample, const size_t pos_ptr, const size_t max_len) -> size_t {
            //verbose("computing MLQ for:  ", pos_sample, " ", pos_pattern);
            #ifdef MEASURE_TIME
            Stopwatch s;
            #endif
            auto lce =  ((pos_sample + 1) >= n) ? 0 : lceToRBounded(slp, pos_sample + 1, pos_ptr, max_len);
            #ifdef MEASURE_TIME
            time_lce += s.seconds();
            ++count_lce_total;
            #endif
            //verbose("computed LCE of:  ", lce, " ", max_len);
            return min(max_len, lce);
        };

        // auto empty_stack = [&] () {
        //     // we only compute (lce_cnt - 1) LCEs as the last one already was computed  
        //     for (int j = 0; j < (lce_cnt - 1); j++) {
        //         if ( (direction[j] == forbidden) ||
        //             ((direction[j] == up) && !thr_lce.skip_preceding_lce(stored_run[j], last_len)) ||
        //             ((direction[j] == down) && !thr_lce.skip_succeeding_lce(stored_run[j], last_len))) {
        //           last_len = compute_lce(stored_sample_pos[j], stored_ptr[j], last_len);
        //         }
        //         write_len_segment(max(1, stored_it[j+1] - stored_it[j]));
        //     }
        // };

        // auto try_skip_lces = [&] (const size_t lces_to_skip) {
        //     const size_t last_delay = lces_to_skip - 1;
        //     const size_t skipped_steps = (stored_it[last_delay] - stored_it[0]);
        //     const int new_r_bound = last_len + skipped_steps; // if we are processing the same MEM, this is its length

        //     const size_t lce = compute_mlq(stored_sample_pos[last_delay], m-stored_it[last_delay]);

        //     lce_is_paused = false;
        //     if (lce == new_r_bound) {
        //         //verbose("first delayed unbounded lce: ", compute_lce(stored_sample_pos[0], m-stored_it[0], 1000));
        //         //verbose("last  delayed unbounded lce: ", compute_lce(stored_sample_pos[last_delay], m-stored_it[last_delay], 1000));
        //         //verbose("lce: ", lce, " last_len+skipped_steps: ", last_len + skipped_steps, " last_len", last_len, " skipped_steps: ", skipped_steps);
        //         //verbose("stored_sample_pos:", stored_sample_pos[last_delay]);
        //         // if we still process the same MEM, we know the lens right ahead
        //         write_len_segment(skipped_steps+1);
        //     } else {
        //         empty_stack();
        //         write_len(lce+1, lce_is_paused);
        //     }
        //     lce_cnt = 0;
        // };

        for (size_t i = 1; i < m; ++i) {
            const auto c = pattern_at(m - i - 1);
            const size_t number_of_runs_of_c = this->bwt.number_of_letter(c);

            if (number_of_runs_of_c == 0) {
                skip_queries = min_ms_len;
            } else if (pos < this->bwt.size() && this->bwt[pos] == c) {
                if (!skip_queries) {
                    write_len(last_len+1);
                    write_ref(last_ref-1);
                }
            } else { // we jump
                const ri::ulint rank = this->bwt.rank(pos, c);

                size_t sa0, sa1; // SA indices of the previous and next entries with BWT[i] = c
                ri::ulint run1, run0;

                if (rank > 0) {
                    sa0 = this->bwt.select(rank-1, c);
                }

                if (rank < number_of_runs_of_c) {
                    sa1 = this->bwt.select(rank, c);
                }

                run1 = this->bwt.run_of_position(sa1);

                const Triplet jump = [&] () -> Triplet {
                    if (rank == 0) {
                        if (skip_queries) {
                            return {sa1, 0, 0};
                        }
                        else if (do_mlq) {
                            const size_t textposStart = subsamples_start(run1, sa1);
                            size_t match_len = compute_mlq(textposStart, m - i - 1);
                            if (match_len >= min_ms_len) {
                                do_mlq = false;
                            }
                            else {
                                skip_queries = min_ms_len - match_len;
                            }
                            return {sa1, textposStart, match_len};
                        }
                        else {
                            const size_t textposStart = subsamples_start(run1, sa1);
                            return {sa1, textposStart, compute_lce(textposStart, last_ref, last_len)};
                        }
                    } else if(rank >= number_of_runs_of_c) {
                        const ri::ulint run0 = this->bwt.run_of_position(sa0);
                        if (skip_queries) {
                            return {sa0, 0, 0};
                        }
                        else if (do_mlq) {
                            const size_t textposLast = subsamples_last(run0, sa0);
                            size_t match_len = compute_mlq(textposLast, m - i - 1);
                            if (match_len >= min_ms_len) {
                                do_mlq = false;
                            }
                            else {
                                skip_queries = min_ms_len - match_len;
                            }
                            return {sa0, textposLast, match_len};
                        }
                        else {
                            const size_t textposLast = subsamples_start(run0, sa0);
                            return {sa0, textposLast, compute_lce(textposLast, last_ref, last_len)};
                        }
                    }
                    // Check thresholds and boundary LCEs first
                    const size_t thr = thresholds[run1];
                    if (pos < thr) {
                        const ri::ulint run0 = this->bwt.run_of_position(sa0);
                        if (!do_mlq && !skip_queries && thr_lce.skip_preceding_lce(run1, last_len)) {
                            const ri::ulint run0 = this->bwt.run_of_position(sa0);
                            const size_t ref0 = subsamples_last(run0, sa0);
                            const size_t len0 = last_len;
                            return {sa0, ref0, len0};
                        } else {
                            if (skip_queries) {
                                return {sa0, 0, 0};
                            }
                            else if (do_mlq) {
                                const size_t textposLast = subsamples_last(run0, sa0);
                                size_t match_len = compute_mlq(textposLast, m - i - 1);
                                if (match_len >= min_ms_len) {
                                    do_mlq = false;
                                }
                                else {
                                    skip_queries = min_ms_len - match_len;
                                }
                                return {sa0, textposLast, match_len};
                            }
                            else {
                                const size_t textposLast = subsamples_start(run0, sa0);
                                return {sa0, textposLast, compute_lce(textposLast, last_ref, last_len)};
                            }
                        }
                    } else {
                        if (!do_mlq && !skip_queries && thr_lce.skip_succeeding_lce(run1, last_len)) {
                            const size_t ref1 = subsamples_start(run1, sa1);
                            const size_t len1 = last_len;
                            return {sa1, ref1, len1};
                        } else {
                            if (skip_queries) {
                                return {sa1, 0, 0};
                            }
                            else if (do_mlq) {
                                const size_t textposStart = subsamples_start(run1, sa1);
                                size_t match_len = compute_mlq(textposStart, m - i - 1);
                                if (match_len >= min_ms_len) {
                                    do_mlq = false;
                                }
                                else {
                                    skip_queries = min_ms_len - match_len;
                                }
                                return {sa1, textposStart, match_len};
                            }
                            else {
                                const size_t textposStart = subsamples_start(run1, sa1);
                                return {sa1, textposStart, compute_lce(textposStart, last_ref, last_len)};
                            }
                        }
                    }
                }();

                if (!skip_queries) {
                    write_len(1 + jump.len);
                    write_ref(jump.ref);
                }
                pos = jump.sa;
            }
            pos = LF(pos, c); //! Perform one backward step
            --skip_queries;
        }
        //verbose("lens_size = ", lens.size(), " refs.size = ", refs.size());
        #ifdef MEASURE_TIME
		cout << "Time backwardsearch: " << time_backwardstep << std::endl;
		cout << "Time lce: " << time_lce << std::endl;
        cout << "Time mlq: " << time_mlq << std::endl;
        cout << "Count mlq: " << count_mlq_total << std::endl;
		cout << "Count lce: " << count_lce_total << std::endl;
		cout << "Count lce skips: " << count_lce_skips << std::endl;
        #endif
        return std::make_pair(lens, refs);
    }

    /*
     * \param i position in the BWT
     * \param c character
     * \return lexicographic rank of cw in bwt
     */
    ulint LF(ri::ulint i, ri::uchar c)
    {
        // //if character does not appear in the text, return empty pair
        // if ((c == 255 and this->F[c] == this->bwt_size()) || this->F[c] >= this->F[c + 1])
        //     return {1, 0};
        //number of c before the interval
        ri::ulint c_before = this->bwt.rank(i, c);
        // number of c inside the interval rn
        ri::ulint l = this->F[c] + c_before;
        return l;
    }

    /**
    * @brief Access the subsampled version of samples_last for testing purposes
    *
    * @param run The run index used for accessing samples_last.
    * @return ulint - The corresponding value of samples_last.
    */
    ulint subsamples_last(const ulint& run) {
        return subsamples_last(run, this->bwt.run_range(run).second);
    }

    /**
    * @brief Access the subsampled version of samples_last
    *
    * @param run The run index used for accessing samples_last.
    * @param pos The last position in the BWT corresponding to the run of interest.
    * @return ulint - The corresponding value of samples_last.
    */
    ulint subsamples_last(const ulint& run, const ulint& pos) {
        DCHECK_EQ(pos, this->bwt.run_range(run).second);

        ON_DEBUG(LF_pos = std::queue<ulint>());
        LF_next = std::queue<ulint>();

        // Check if the run is subsampled using subsampled_last_samples_bv
        if (subsampled_last_samples_bv[run] == 1) {
            return this->samples_last[subsampled_last_samples_bv.rank(run)];
        }

        auto n = this->bwt_size();

        // Find the character at current_pos
        auto c = this->bwt[pos];

        ON_DEBUG(LF_pos.push(pos));

        // Traverse the BWT to find a subsampled run
        ulint current_pos = this->LF(pos, c);
        ulint current_run = this->bwt.run_of_position(current_pos);
        size_t k = 1;

        LF_next.push(current_pos);

        while (subsampled_last_samples_bv[current_run] == 0 ||
            (current_pos != n - 1 && current_run == this->bwt.run_of_position(current_pos + 1))) { 
            
            ON_DEBUG(LF_pos.push(current_pos));

            c = this->bwt[current_pos];
            current_pos = this->LF(current_pos, c);
            current_run = this->bwt.run_of_position(current_pos);
            k++;

            LF_next.push(current_pos);
        }

        assert(k < this->maxLF);

        return this->samples_last[subsampled_last_samples_bv.rank(current_run)] + k;
    }


    /**
    * @brief Access the subsampled version of samples_start for testing purposes
    *
    * @param run The run index used for accessing samples_start.
    * @return ulint - The corresponding value of samples_start.
    */
    ulint subsamples_start(const ulint& run) {
        return subsamples_start(run, this->bwt.run_range(run).first);
    }


    /**
    * @brief Access the subsampled version of samples_start
    *
    * @param run The run index used for accessing samples_start.
    * @param pos The first position in the BWT corresponding to the run of interest.
    * @return ulint - The corresponding value of samples_start.
    */
    ulint subsamples_start(const ulint& run, const ulint& pos) {
        DCHECK_EQ(pos, this->bwt.run_range(run).first);

        ON_DEBUG(LF_pos = std::queue<ulint>());
        LF_next = std::queue<ulint>();

        // Check if the run is subsampled using subsampled_start_samples_bv
        if (subsampled_start_samples_bv[run] == 1) {
            return this->samples_start[subsampled_start_samples_bv.rank(run)];
        }

        // Find the character at current_pos
        auto c = this->bwt[pos];

        ON_DEBUG(LF_pos.push(pos));

        // Traverse the BWT to find a subsampled run
        ulint current_pos = this->LF(pos, c);
        ulint current_run = this->bwt.run_of_position(current_pos);
        size_t k = 1;

        LF_next.push(current_pos);

        while (subsampled_start_samples_bv[current_run] == 0 ||
            (current_pos != 0 && current_run == this->bwt.run_of_position(current_pos - 1))) { 
            
            ON_DEBUG(LF_pos.push(current_pos));

            c = this->bwt[current_pos];
            current_pos = this->LF(current_pos, c);
            current_run = this->bwt.run_of_position(current_pos);
            k++;

            LF_next.push(current_pos);
        }

        assert(k < this->maxLF);

        return this->samples_start[subsampled_start_samples_bv.rank(current_run)] + k;
    }

    /* serialize the structure to the ostream
     * \param out     the ostream
     */
    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") // const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;

        out.write((char *)&this->terminator_position, sizeof(this->terminator_position));
        written_bytes += sizeof(this->terminator_position);
        written_bytes += my_serialize(this->F, out, child, "F");
        written_bytes += this->bwt.serialize(out);
        written_bytes += this->samples_last.serialize(out);

        written_bytes += subsampled_last_samples_bv.serialize(out);

        written_bytes += samples_start.serialize(out, child, "samples_start");

        written_bytes += subsampled_start_samples_bv.serialize(out);

        written_bytes += thresholds.serialize(out, child, "thresholds");
        written_bytes += thr_lce.serialize(out, child, "thr_lce");

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    /* load the structure from the istream
     * \param in the istream
     */
    void load(std::istream &in, const std::string& filename, const size_t& _maxLF)
    {
        in.read((char *)&this->terminator_position, sizeof(this->terminator_position));
        my_load(this->F, in);
        this->bwt.load(in);
        this->r = this->bwt.number_of_runs();
        this->samples_last.load(in);
        subsampled_last_samples_bv.load(in);
        this->samples_start.load(in);
        subsampled_start_samples_bv.load(in);

        thresholds.load(in, &this->bwt);
        thr_lce.load(in, &this->bwt);
        
        load_grammar(filename);
    }

    // // From r-index
    // ulint get_last_run_sample()
    // {
    //     return (samples_last[r - 1] + 1) % bwt.size();
    // }

protected :

    vector<ulint> build_F_(std::ifstream &heads, std::ifstream &lengths)
    {
        heads.clear();
        heads.seekg(0);
        lengths.clear();
        lengths.seekg(0);

        this->F = vector<ulint>(256, 0);
        int c;
        {
            ulint i = 0;
            while ((c = heads.get()) != EOF) {
                size_t length;
                lengths.read((char *)&length, 5);
                if (c > TERMINATOR)
                    this->F[c] += length;
                else {
                    this->F[TERMINATOR] += length;
                    this->terminator_position = i;
                }
                i++;
            }
        }
        for (ulint i = 255; i > 0; --i)
            this->F[i] = this->F[i - 1];
        this->F[0] = 0;
        for (ulint i = 1; i < 256; ++i)
            this->F[i] += this->F[i - 1];
        return this->F;
    }
};

#endif /* end of include guard: _MS_POINTERS_HH */
