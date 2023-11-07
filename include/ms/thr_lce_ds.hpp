/* thr_lce_ds - Stores threshold LCE values in compressed and plain ways 
    Copyright (C) 2022 Nathaniel Brown

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
   \file thr_lce_ds.hpp
   \brief thr_lce_ds.hpp Stores LCE values between run boundaries and thresholds in compressed and plain ways.
   \author Nathaniel Brown
   \date 10/28/2022
*/

#ifndef _MS_THR_LCE_DS_HH
#define _MS_THR_LCE_DS_HH

#include <common.hpp>
#include <cmath>

//#include <malloc_count.h>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/dac_vector.hpp>

#include <ms_rle_string.hpp>

// Used to represent the 'stored' vectors
template < class bv_t = sd_vector<> >
class rank_bv
{
public:
    typedef typename bv_t::rank_1_type bv_rank;

    bv_t bits;
    bv_rank bit_rank;

    rank_bv() {}

    rank_bv(vector<bool> vec) {
        // Convert boolean vector to specified bit vector
        sdsl::bit_vector bv(vec.size());

        for (size_t i = 0; i < vec.size(); ++i)
            bv[i] = vec[i];

        bits =  bv_t(bv);
        bit_rank = bv_rank(&bits);
    }

    size_t operator[] (size_t& i)
    {
        return bits[i];
    }

    size_t rank(size_t i) 
    {
        return bit_rank(i);
    }

    /* serialize the structure to the ostream
    * \param out     the ostream
    */
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;

        written_bytes += bits.serialize(out, v, "rank_bv");

        return written_bytes;
    }

    /* load the structure from the istream
    * \param in the istream
    */
    void load(std::istream &in)
    {
        bits.load(in);
        bit_rank = bv_rank(&bits);
    }
};

// TODO: Add template/constructor for DAC on top of int_vector (a vector type), might need a wrapper
template <class rle_string_t = ms_rle_string_sd>
class thr_lce_plain
{
public:
    int max_bits; // Width of maximum bits we store, 0 to use as large as necessary (lg n)
    ulint max_size; // maximum unsigned value given by an integer of max_bits width

    int_vector<> before_thr_lce;
    int_vector<> after_thr_lce;

    rle_string_t *bwt;

    typedef size_t size_type;

    thr_lce_plain()
    {
        bwt=nullptr;
        max_bits=0;
    }

    thr_lce_plain(std::string filename, rle_string_t* bwt_, int byte_width = 0):bwt(bwt_)
    {
        int bit_width = 8*byte_width;
        int log_n = bitsize(uint64_t(bwt->size()));
        int vec_width = log_n;
        max_bits = 0;

        if (bit_width && bit_width < vec_width)
        {
            max_bits = bit_width;
            vec_width = max_bits;
            max_size = pow(2, max_bits) - 1;
        }

        verbose("Reading threshold/boundary LCE values from file");

        std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

        std::string tmp_filename = filename + std::string(".thr_lce");

        struct stat filestat;
        FILE *fd;

        if ((fd = fopen(tmp_filename.c_str(), "r")) == nullptr)
            error("open() file " + tmp_filename + " failed");

        int fn = fileno(fd);
        if (fstat(fn, &filestat) < 0)
            error("stat() file " + tmp_filename + " failed");

        if (filestat.st_size % THRBYTES != 0)
            error("invalid file " + tmp_filename);

        assert(filestat.st_size / THRBYTES == 2*bwt->number_of_runs());
        size_t length = bwt->number_of_runs();

        before_thr_lce = int_vector<>(length, 0, vec_width);
        after_thr_lce = int_vector<>(length, 0, vec_width);

        for (size_t i = 0; i < length; ++i)
        {
            size_t before_lce = 0;
            size_t after_lce = 0;

            if ((fread(&before_lce, THRBYTES, 1, fd)) != 1)
                error("fread() file " + tmp_filename + " failed");
            if ((fread(&after_lce, THRBYTES, 1, fd)) != 1)
                error("fread() file " + tmp_filename + " failed");

            // If using capped bit size, store the max to signal an overflow
            if (max_bits) 
            {
                before_thr_lce[i] = std::max<uint64_t>(before_lce, max_size);
                after_thr_lce[i] = std::max<uint64_t>(after_lce, max_size);
            }
            else
            {
                before_thr_lce[i] = before_lce;
                after_thr_lce[i] = after_lce;
            }
        }

        fclose(fd);

        std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

        //verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
    }

    // Destructor
    ~thr_lce_plain() 
    {
       // NtD
    }

    // Copy constructor
    thr_lce_plain(const thr_lce_plain &other)
        :before_thr_lce(other.before_thr_lce),
        after_thr_lce(other.after_thr_lce),
        max_bits(other.max_bits),
        bwt(other.bwt)
    {
    }

    friend void swap(thr_lce_plain &first, thr_lce_plain &second) // nothrow
    {
        using std::swap;

        swap(first.before_thr_lce, second.before_thr_lce);
        swap(first.after_thr_lce, second.after_thr_lce);
        swap(first.max_bits, second.max_bits);
        swap(first.bwt, second.bwt);
    }

    // Copy assignment
    thr_lce_plain &operator=(thr_lce_plain other) 
    {
        swap(*this,other);
        
        return *this;
    }

    // Move constructor
    thr_lce_plain(thr_lce_plain &&other) noexcept
        : thr_lce_plain()
    {
        swap(*this, other);
    }

    // std::pair<ulint, ulint> operator[] (size_t& i)
    // {
    //     assert( i < before_thr_lce.size());
    //     assert( i < after_thr_lce.size());
    //     return std::pair<ulint, ulint>(before_thr_lce[i], after_thr_lce[i]);
    // }

    bool skip_preceding_lce(const size_t run, const size_t length)
    {
        // If using capped bit size, skip if we overflow
        if (max_bits)
        {
            const size_t thr_lce = before_thr_lce[run];
            return thr_lce < max_size && length <= thr_lce;
        }

        return length <= before_thr_lce[run];
    }

    bool skip_succeeding_lce(const size_t run, const size_t length)
    {
        // If using capped bit size, skip if we overflow
        if (max_bits)
        {
            uint64_t thr_lce = after_thr_lce[run];
            return thr_lce < max_size && length <= thr_lce;
        }

        return length <= after_thr_lce[run];
    }

    /* serialize the structure to the ostream
     * \param out     the ostream
     */
    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") // const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;

        out.write((char *)&max_bits, sizeof(max_bits));
        written_bytes += sizeof(max_bits);

        written_bytes += before_thr_lce.serialize(out, child, "before_thr_lce");
        written_bytes += after_thr_lce.serialize(out, child, "after_thr_lce");

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    /* load the structure from the istream
     * \param in the istream
     */
    void load(std::istream &in, rle_string_t *bwt_)
    {
        in.read((char *)&max_bits, sizeof(max_bits));
        max_size = (max_bits) ? pow(2, max_bits) - 1 : 0;

        before_thr_lce.load(in);
        after_thr_lce.load(in);

        bwt = bwt_;
    }

    std::string get_file_extension() const
    {
        return ".thr_lce_plain_ds";
    }
};

// TODO: Add template for vector
template <class rle_string_t = ms_rle_string_sd,
          class bv_t = sd_vector<> >
class thr_lce_bv
{
public:
    rank_bv<bv_t> preceding_stored; // we mark non zero positions, and those below overflow
    rank_bv<bv_t> succeeding_stored; // we mark non zero positions, and those below overflow

    int max_bits; // Width of maximum bits we store, 0 to use as large as necessary (lg n)
    ulint max_size; // maximum unsigned value given by an integer of max_bits width

    int_vector<> before_thr_lce;
    int_vector<> after_thr_lce;

    rle_string_t *bwt;

    typedef size_t size_type;

    thr_lce_bv()
    {
        bwt=nullptr;
        max_bits=0;
    }

    thr_lce_bv(std::string filename, rle_string_t* bwt_, int bit_width = 0):bwt(bwt_)
    {
        int log_n = bitsize(uint64_t(bwt->size()));
        int vec_width = log_n;
        max_bits = 0;

        if (bit_width && bit_width < vec_width)
        {
            max_bits = bit_width;
            vec_width = max_bits;
            max_size = max_bits^2 - 1;
        }

        verbose("Reading threshold/boundary LCE values from file");

        std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

        std::vector<ulint> stored_thr_lce_before = std::vector<ulint>();
        std::vector<ulint> stored_thr_lce_after = std::vector<ulint>();

        std::string tmp_filename = filename + std::string(".thr_lce");

        struct stat filestat;
        FILE *fd;

        if ((fd = fopen(tmp_filename.c_str(), "r")) == nullptr)
            error("open() file " + tmp_filename + " failed");

        int fn = fileno(fd);
        if (fstat(fn, &filestat) < 0)
            error("stat() file " + tmp_filename + " failed");

        if (filestat.st_size % THRBYTES != 0)
            error("invalid file " + tmp_filename);

        assert(filestat.st_size / THRBYTES == 2*bwt->number_of_runs());
        size_t length = bwt->number_of_runs();

        std::vector<bool> stored_positions_before(length, false);
        std::vector<bool> stored_positions_after(length, false);

        for (size_t i = 0; i < length; ++i)
        {
            size_t before_lce = 0;
            size_t after_lce = 0;

            if ((fread(&before_lce, THRBYTES, 1, fd)) != 1)
                error("fread() file " + tmp_filename + " failed");
            if ((fread(&after_lce, THRBYTES, 1, fd)) != 1)
                error("fread() file " + tmp_filename + " failed");

            // If using capped bit size, don't mark the bit and store, or if the value is 0
            // Also, logic is convoluted, should rewrite (TODO)
            if (!max_bits || before_lce < max_size) 
            {
                if (before_lce > 0)
                {
                    stored_thr_lce_before.push_back(before_lce);
                    stored_positions_before[i] = true;
                }
            }

            if (!max_bits || after_lce < max_size) 
            {
                if (after_lce > 0)
                {
                    stored_thr_lce_after.push_back(after_lce);
                    stored_positions_after[i] = true;
                }
            }
        }

        fclose(fd);

        // Do a naive copy to make sure the vec_width is correct (better way?)
        before_thr_lce = int_vector<>(stored_thr_lce_before.size(), 0, vec_width);
        preceding_stored = rank_bv<>(stored_positions_before);
        assert(preceding_stored.rank(stored_positions_before.size()) == stored_thr_lce_before.size());

        for (size_t i = 0; i < stored_thr_lce_before.size(); ++i)
        {
            before_thr_lce[i] = stored_thr_lce_before[i];
        }

        // Do a naive copy to make sure the vec_width is correct (better way?)
        after_thr_lce = int_vector<>(stored_thr_lce_after.size(), 0, vec_width);
        succeeding_stored = rank_bv<>(stored_positions_after);
        assert(succeeding_stored.rank(stored_positions_after.size()) == stored_thr_lce_after.size());

        for (size_t i = 0; i < stored_thr_lce_after.size(); ++i)
        {
            after_thr_lce[i] = stored_thr_lce_after[i];
        }

        std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

        //verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
    }

    // Destructor
    ~thr_lce_bv() 
    {
       // NtD
    }

    // Copy constructor
    thr_lce_bv(const thr_lce_bv &other)
        :before_thr_lce(other.before_thr_lce),
        after_thr_lce(other.after_thr_lce),
        preceding_stored(other.preceding_stored),
        succeeding_stored(other.succeeding_stored),
        max_bits(other.max_bits),
        bwt(other.bwt)
    {
    }

    friend void swap(thr_lce_bv &first, thr_lce_bv &second) // nothrow
    {
        using std::swap;

        swap(first.before_thr_lce, second.before_thr_lce);
        swap(first.after_thr_lce, second.after_thr_lce);
        swap(first.max_bits, second.max_bits);
        swap(first.preceding_stored, second.preceding_stored);
        swap(first.succeeding_stored, second.succeeding_stored);
        swap(first.bwt, second.bwt);
    }

    // Copy assignment
    thr_lce_bv &operator=(thr_lce_bv other) 
    {
        swap(*this,other);
        
        return *this;
    }

    // Move constructor
    thr_lce_bv(thr_lce_bv &&other) noexcept
        : thr_lce_bv()
    {
        swap(*this, other);
    }

    // std::pair<ulint, ulint> operator[] (size_t& i)
    // {
    //     assert( i < before_thr_lce.size());
    //     assert( i < after_thr_lce.size());
    //     return std::pair<ulint, ulint>(before_thr_lce[i], after_thr_lce[i]);
    // }

    bool skip_preceding_lce(size_t run, size_t length)
    {
        if (!preceding_stored[run]) return false;

        size_t lce_idx = preceding_stored.rank(run);
        return length <= before_thr_lce[lce_idx];
    }

    bool skip_succeeding_lce(size_t run, size_t length)
    {
        if (!succeeding_stored[run]) return false;

        size_t lce_idx = succeeding_stored.rank(run);
        return length <= after_thr_lce[lce_idx];
    }

    /* serialize the structure to the ostream
     * \param out     the ostream
     */
    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") // const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;

        out.write((char *)&max_bits, sizeof(max_bits));
        written_bytes += sizeof(max_bits);

        written_bytes += before_thr_lce.serialize(out, child, "before_thr_lce");
        written_bytes += after_thr_lce.serialize(out, child, "after_thr_lce");

        written_bytes += preceding_stored.serialize(out, child, "preceding_stored");
        written_bytes += succeeding_stored.serialize(out, child, "succeeding_stored");

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    /* load the structure from the istream
     * \param in the istream
     */
    void load(std::istream &in, rle_string_t *bwt_)
    {
        in.read((char *)&max_bits, sizeof(max_bits));
        max_size = (max_bits) ? max_bits^2 - 1 : 0;

        before_thr_lce.load(in);
        after_thr_lce.load(in);
        
        preceding_stored.load(in);
        succeeding_stored.load(in);
        bwt = bwt_;
    }

    std::string get_file_extension() const
    {
        return ".thr_lce_bv_ds";
    }
};

template <class rle_string_t = ms_rle_string_sd,
          class vec_t = dac_vector_dp<> >
class thr_lce_dac
{
public:
    vec_t before_thr_lce;
    vec_t after_thr_lce;

    rle_string_t *bwt;

    typedef size_t size_type;

    thr_lce_dac()
    {
        bwt=nullptr;
    }

    thr_lce_dac(std::string filename, rle_string_t* bwt_, int bytes = 0):bwt(bwt_)
    {
        verbose("Reading threshold/boundary LCE values from file");

        std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

        std::string tmp_filename = filename + std::string(".thr_lce");

        struct stat filestat;
        FILE *fd;

        if ((fd = fopen(tmp_filename.c_str(), "r")) == nullptr)
            error("open() file " + tmp_filename + " failed");

        int fn = fileno(fd);
        if (fstat(fn, &filestat) < 0)
            error("stat() file " + tmp_filename + " failed");

        if (filestat.st_size % THRBYTES != 0)
            error("invalid file " + tmp_filename);

        assert(filestat.st_size / THRBYTES == 2*bwt->number_of_runs());
        size_t length = bwt->number_of_runs();

        std::vector<ulint> before_vals = std::vector<ulint>(length, 0);
        std::vector<ulint> after_vals = std::vector<ulint>(length, 0);

        for (size_t i = 0; i < length; ++i)
        {
            size_t before_lce = 0;
            size_t after_lce = 0;

            if ((fread(&before_lce, THRBYTES, 1, fd)) != 1)
                error("fread() file " + tmp_filename + " failed");
            if ((fread(&after_lce, THRBYTES, 1, fd)) != 1)
                error("fread() file " + tmp_filename + " failed");

            before_vals[i] = before_lce;
            after_vals[i] = after_lce;
        }

        before_thr_lce = vec_t(before_vals);
        after_thr_lce = vec_t(after_vals);

        fclose(fd);

        std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

        //verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
    }

    // Destructor
    ~thr_lce_dac() 
    {
       // NtD
    }

    // Copy constructor
    thr_lce_dac(const thr_lce_dac &other)
        :before_thr_lce(other.before_thr_lce),
        after_thr_lce(other.after_thr_lce),
        bwt(other.bwt)
    {
    }

    friend void swap(thr_lce_dac &first, thr_lce_dac &second) // nothrow
    {
        using std::swap;

        swap(first.before_thr_lce, second.before_thr_lce);
        swap(first.after_thr_lce, second.after_thr_lce);
        swap(first.bwt, second.bwt);
    }

    // Copy assignment
    thr_lce_dac &operator=(thr_lce_dac other) 
    {
        swap(*this,other);
        
        return *this;
    }

    // Move constructor
    thr_lce_dac(thr_lce_dac &&other) noexcept
        : thr_lce_dac()
    {
        swap(*this, other);
    }

    bool skip_preceding_lce(const size_t run, const size_t length)
    {
        return length <= before_thr_lce[run];
    }

    bool skip_succeeding_lce(const size_t run, const size_t length)
    {
        return length <= after_thr_lce[run];
    }

    /* serialize the structure to the ostream
     * \param out     the ostream
     */
    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") // const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;

        written_bytes += before_thr_lce.serialize(out, child, "before_thr_lce");
        written_bytes += after_thr_lce.serialize(out, child, "after_thr_lce");

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    /* load the structure from the istream
     * \param in the istream
     */
    void load(std::istream &in, rle_string_t *bwt_)
    {
        before_thr_lce.load(in);
        after_thr_lce.load(in);

        bwt = bwt_;
    }

    std::string get_file_extension() const
    {
        return ".thr_lce_dac_ds";
    }
};

// TODO: Add template for vector
template <class rle_string_t = ms_rle_string_sd,
          class vec_t = dac_vector_dp<>, 
          class bv_t = sd_vector<> >
class thr_lce_bv_dac
{
public:
    rank_bv<bv_t> preceding_stored; // we mark non zero positions, and those below overflow
    rank_bv<bv_t> succeeding_stored; // we mark non zero positions, and those below overflow

    vec_t before_thr_lce;
    vec_t after_thr_lce;

    rle_string_t *bwt;

    typedef size_t size_type;

    thr_lce_bv_dac()
    {
        bwt=nullptr;
    }

    thr_lce_bv_dac(std::string filename, rle_string_t* bwt_, int bytes = 0):bwt(bwt_)
    {
        verbose("Reading threshold/boundary LCE values from file");

        std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

        std::vector<ulint> stored_thr_lce_before = std::vector<ulint>();
        std::vector<ulint> stored_thr_lce_after = std::vector<ulint>();

        std::string tmp_filename = filename + std::string(".thr_lce");

        struct stat filestat;
        FILE *fd;

        if ((fd = fopen(tmp_filename.c_str(), "r")) == nullptr)
            error("open() file " + tmp_filename + " failed");

        int fn = fileno(fd);
        if (fstat(fn, &filestat) < 0)
            error("stat() file " + tmp_filename + " failed");

        if (filestat.st_size % THRBYTES != 0)
            error("invalid file " + tmp_filename);

        assert(filestat.st_size / THRBYTES == 2*bwt->number_of_runs());
        size_t length = bwt->number_of_runs();

        std::vector<bool> stored_positions_before(length, false);
        std::vector<bool> stored_positions_after(length, false);

        for (size_t i = 0; i < length; ++i)
        {
            size_t before_lce = 0;
            size_t after_lce = 0;

            if ((fread(&before_lce, THRBYTES, 1, fd)) != 1)
                error("fread() file " + tmp_filename + " failed");
            if ((fread(&after_lce, THRBYTES, 1, fd)) != 1)
                error("fread() file " + tmp_filename + " failed");

            if (before_lce > 0)
            {
                stored_thr_lce_before.push_back(before_lce);
                stored_positions_before[i] = true;
            }

            if (after_lce > 0)
            {
                stored_thr_lce_after.push_back(after_lce);
                stored_positions_after[i] = true;
            }
        }

        fclose(fd);

        before_thr_lce = vec_t(stored_thr_lce_before);
        preceding_stored = rank_bv<>(stored_positions_before);
        assert(preceding_stored.rank(stored_positions_before.size()) == stored_thr_lce_before.size());

        after_thr_lce = vec_t(stored_thr_lce_after);
        succeeding_stored = rank_bv<>(stored_positions_after);
        assert(succeeding_stored.rank(stored_positions_after.size()) == stored_thr_lce_after.size());

        std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

        //verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
    }

    // Destructor
    ~thr_lce_bv_dac() 
    {
       // NtD
    }

    // Copy constructor
    thr_lce_bv_dac(const thr_lce_bv_dac &other)
        :before_thr_lce(other.before_thr_lce),
        after_thr_lce(other.after_thr_lce),
        preceding_stored(other.preceding_stored),
        succeeding_stored(other.succeeding_stored),
        bwt(other.bwt)
    {
    }

    friend void swap(thr_lce_bv_dac &first, thr_lce_bv_dac &second) // nothrow
    {
        using std::swap;

        swap(first.before_thr_lce, second.before_thr_lce);
        swap(first.after_thr_lce, second.after_thr_lce);
        swap(first.preceding_stored, second.preceding_stored);
        swap(first.succeeding_stored, second.succeeding_stored);
        swap(first.bwt, second.bwt);
    }

    // Copy assignment
    thr_lce_bv_dac &operator=(thr_lce_bv_dac other) 
    {
        swap(*this,other);
        
        return *this;
    }

    // Move constructor
    thr_lce_bv_dac(thr_lce_bv_dac &&other) noexcept
        : thr_lce_bv_dac()
    {
        swap(*this, other);
    }

    bool skip_preceding_lce(size_t run, size_t length)
    {
        if (!preceding_stored[run]) return false;

        size_t lce_idx = preceding_stored.rank(run);
        return length <= before_thr_lce[lce_idx];
    }

    bool skip_succeeding_lce(size_t run, size_t length)
    {
        if (!succeeding_stored[run]) return false;

        size_t lce_idx = succeeding_stored.rank(run);
        return length <= after_thr_lce[lce_idx];
    }

    /* serialize the structure to the ostream
     * \param out     the ostream
     */
    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") // const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;

        written_bytes += before_thr_lce.serialize(out, child, "before_thr_lce");
        written_bytes += after_thr_lce.serialize(out, child, "after_thr_lce");

        written_bytes += preceding_stored.serialize(out, child, "preceding_stored");
        written_bytes += succeeding_stored.serialize(out, child, "succeeding_stored");

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    /* load the structure from the istream
     * \param in the istream
     */
    void load(std::istream &in, rle_string_t *bwt_)
    {
        before_thr_lce.load(in);
        after_thr_lce.load(in);
        
        preceding_stored.load(in);
        succeeding_stored.load(in);
        bwt = bwt_;
    }

    std::string get_file_extension() const
    {
        return ".thr_lce_bv_dac_ds";
    }
};

#endif /* end of include guard: _MS_THR_LCE_DS_HH */
