#include <seqan3/argument_parser/all.hpp>     // includes all necessary headers
#include <seqan3/io/alignment_file/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/range/views/all.hpp>
// #include <seqan3/range/views/all.hpp>
#include <seqan3/core/debug_stream.hpp>  // our custom output stream
#include <seqan3/std/filesystem>          // use std::filesystem::path
#include <seqan3/std/ranges>
#include <time.h>

using namespace seqan3;

struct reference_storage_t
{
    std::vector<std::string> ids;
    std::vector<std::vector<dna5>> seqs;
};

void read_reference(std::filesystem::path const & reference_path,
                    reference_storage_t & ref_storage)
{
    debug_stream << "INFO: reading reference file" << std::endl;
    time_t start_time = time(NULL);

    sequence_file_input reference_in{reference_path, fields<field::ID, field::SEQ>{}};
    for (auto & [id, seq] : reference_in)
    {
        ref_storage.ids.push_back(std::move(id));
        ref_storage.seqs.push_back(std::move(seq));
    }

    debug_stream << "INFO: finished reading reference file (" << ref_storage.seqs.size() << " records) " 
                 << "in " << time(NULL)-start_time << " seconds" << std::endl;
}

struct single_ref_align_storage_t
{
    std::int32_t ref_idx{-1};
    std::vector<std::string> read_ids;
    std::vector<std::vector<dna5>> read_seqs;
    std::vector<std::uint16_t> flags;
    std::vector<std::vector<cigar>> align_cigars;
    std::vector<std::int32_t> ref_offsets;
    std::vector<std::int32_t> read_offsets;
};

void print_debug(std::string const & read_id,
                 std::uint16_t const flag,
                 std::optional<int32_t> const & ref_idx,
                 reference_storage_t const & ref_storage,
                 std::optional<int32_t> const & ref_offset,
                 std::vector<cigar> const & align_cigar,
                 std::vector<dna5> const & read_seq,
                 std::int32_t const read_offset)
{
    debug_stream << std::endl;
    debug_stream << "read id:     " << read_id << std::endl;
    debug_stream << "flag:        " << flag << std::endl;
    debug_stream << "is revcomp:  " << ((flag & 0x10) == 0x10) << std::endl;
    debug_stream << "ref idx:     " << ref_idx << std::endl;
    debug_stream << "ref id:      ";
    if (ref_idx) debug_stream << ref_storage.ids[ref_idx.value()];
    else debug_stream << "no ref";
    debug_stream << std::endl;
    debug_stream << "ref offset:  " << ref_offset << std::endl;
    debug_stream << "cigar:       " << align_cigar << std::endl;
    debug_stream << "read seq:    " << read_seq << std::endl;
    // debug_stream << "revcomp seq: " << (read_seq | views::complement | std::views::reverse | views::to<std::vector>) << std::endl;
    debug_stream << "read offset: " << read_offset << std::endl;
    debug_stream << "#################" << std::endl;
}

std::int32_t get_last_base_ref_pos(std::int32_t const begin_ref_offset, 
                                   std::vector<cigar> const & read_cigar)
{
    std::int32_t end_ref_offset{begin_ref_offset};

    for (auto & cigar_op : read_cigar)
    {
        if (cigar_op == 'M'_cigar_op or cigar_op == 'D'_cigar_op) end_ref_offset += get<0>(cigar_op);
    }

    return end_ref_offset - 1;
}

std::int32_t find_alignment_begin_position_on_read(std::vector<cigar> const & read_cigar)
{
        // Checks if the given index in the cigar vector is a soft clip.
        auto soft_clipping_at = [&] (size_t const index) { return read_cigar[index] == 'S'_cigar_op; };
        // Checks if the given index in the cigar vector is a hard clip.
        auto hard_clipping_at = [&] (size_t const index) { return read_cigar[index] == 'H'_cigar_op; };
        // Checks if the given cigar vector as at least min_size many elements.
        auto vector_size_at_least = [&] (size_t const min_size) { return read_cigar.size() >= min_size; };
        // Returns the cigar count of the ith cigar element in the given cigar vector.
        auto cigar_count_at = [&] (size_t const index) { return get<0>(read_cigar[index]); };

        // check for soft clipping at the first two positions
        if (vector_size_at_least(1) && soft_clipping_at(0))
            return cigar_count_at(0);
        else if (vector_size_at_least(2) && hard_clipping_at(0) && soft_clipping_at(1))
            return cigar_count_at(1);
}

std::int32_t find_overlap_begin_position_on_read(std::int32_t const read_i_ref_offset,
                                                 std::int32_t const read_j_ref_offset,
                                                 std::vector<cigar> const & read_i_cigar)
{
    std::int32_t read_i_current_pos_on_read{0};
    std::int32_t read_i_current_pos_on_ref{read_i_ref_offset};

    for (auto & cigar_op : read_i_cigar)
    {
        std::int32_t operation_count = get<0>(cigar_op);

        if (cigar_op == 'S'_cigar_op or cigar_op == 'I'_cigar_op)
        {
            read_i_current_pos_on_read += operation_count;
        }
        else if (cigar_op == 'D'_cigar_op) // DANGER ZONE
        {
            read_i_current_pos_on_ref += operation_count;

            if (read_i_current_pos_on_ref >= read_j_ref_offset) break;
        }
        else if (cigar_op == 'M'_cigar_op)
        {
            std::int32_t offset = std::min(operation_count, read_j_ref_offset - read_i_current_pos_on_ref);
            read_i_current_pos_on_read += offset;
            read_i_current_pos_on_ref += operation_count;
        }

        if (read_i_current_pos_on_ref >= read_j_ref_offset) break;
    }

    return read_i_current_pos_on_read;
}

std::tuple<int32_t, int32_t> find_anchor(std::vector<cigar> const & read_i_cigar,
                                         std::int32_t const read_i_ref_offset,
                                         std::vector<cigar> const & read_j_cigar,
                                         std::int32_t const read_j_ref_offset)
{
    std::int32_t read_j_alignment_begin_position_on_read = find_alignment_begin_position_on_read(read_j_cigar);
    std::int32_t read_i_overlap_begin_position_on_read = find_overlap_begin_position_on_read(read_i_ref_offset, read_j_ref_offset, read_i_cigar);

    return std::make_tuple(read_i_overlap_begin_position_on_read, read_j_alignment_begin_position_on_read);
}

bool perfect_overlap(std::vector<dna5> const & read_i_seq, 
                     std::vector<dna5> const & read_j_seq,
                     std::tuple<int32_t, int32_t> const anchor_positions,
                     std::uint16_t const min_overlap_length)
{
    auto & read_i_anchor_pos = std::get<0>(anchor_positions);
    auto & read_j_anchor_pos = std::get<1>(anchor_positions);

    std::int32_t upstream_overlap_length = std::min(read_i_anchor_pos, read_j_anchor_pos);
    std::int32_t downstream_overlap_length = std::min(read_i_seq.size() - read_i_anchor_pos,
                                                      read_j_seq.size() - read_j_anchor_pos);

    if (upstream_overlap_length + downstream_overlap_length < min_overlap_length) return false;

    auto read_i_overlap_subseq = read_i_seq | views::slice(read_i_anchor_pos - upstream_overlap_length, read_i_anchor_pos + downstream_overlap_length);
    auto read_j_overlap_subseq = read_j_seq | views::slice(read_j_anchor_pos - upstream_overlap_length, read_j_anchor_pos + downstream_overlap_length);

    // debug_stream << read_i_overlap_subseq << std::endl;
    // debug_stream << read_j_overlap_subseq << std::endl;
    // debug_stream << std::endl;

    return std::ranges::equal(read_i_overlap_subseq, read_j_overlap_subseq);
}

bool are_compatible_reads(std::vector<dna5> const & read_i_seq, 
                          std::vector<dna5> const & read_j_seq,
                          std::tuple<int32_t, int32_t> const anchor_positions,
                          std::uint16_t const min_overlap_length)
{
    // Test for perfect identity of the overlap between the 2 reads
    if (perfect_overlap(read_i_seq, read_j_seq, anchor_positions, min_overlap_length)) return true;

    // If we also have the end pos of the alignment on the reads we can also check
    // for length difference to quickly find non compatible reads, even when allowing errors


    // Additionally, we can perform 2 anchored banded alignment (up and downstream the anchor)
    // to check for overlaps with errors


    return false;
}

void write_compatible_reads(std::string const & read_i_id,
                            std::string const & read_j_id,
                            std::filesystem::path const & output_basepath)
{
    std::cout << read_i_id << "\t" << read_j_id << std::endl;
}

void find_reads_overlaps(single_ref_align_storage_t const & align_storage,
                         std::filesystem::path const & output_basepath, 
                         std::uint16_t const min_overlap_length)
{
    // debug_stream << std::endl << std::endl;
    // debug_stream << "Ref idx:  " << align_storage.ref_idx << std::endl;
    // debug_stream << "Read names:  " << std::endl;
    // debug_stream << align_storage.read_ids << std::endl;
    // debug_stream << "#################" << std::endl << std::endl;

    std::size_t max_comparison_nb = 1000000000;

    std::size_t reads_number = align_storage.read_ids.size();

    if(reads_number>1)
    for(size_t i = 0; i < reads_number-1; ++i)
    {
        auto & read_i_id = align_storage.read_ids[i];
        auto & read_i_seq = align_storage.read_seqs[i];
        auto & read_i_cigar = align_storage.align_cigars[i];
        auto & read_i_ref_offset = align_storage.ref_offsets[i];

        std::int32_t read_i_last_base_ref_pos = get_last_base_ref_pos(read_i_ref_offset, read_i_cigar);

        // debug_stream << std::endl;
        // debug_stream << "ref offset:  " << read_i_ref_offset << std::endl;
        // debug_stream << "cigar:       " << read_i_cigar << std::endl;
        // debug_stream << "last base:   " << read_i_last_base_ref_pos << std::endl;

        std::int32_t forward_comparison_count{0};

        for(size_t j = i+1; j < reads_number; ++j)
        {
            auto & read_j_id = align_storage.read_ids[j];
            auto & read_j_seq = align_storage.read_seqs[j];
            auto & read_j_cigar = align_storage.align_cigars[j];
            auto & read_j_ref_offset = align_storage.ref_offsets[j];

            // Stop comparisons when there is no overlap left possible
            if (read_j_ref_offset > read_i_last_base_ref_pos) break;

            // Find anchor positions on reads when we know it exists an overlap of size >= 1
            std::tuple<int32_t, int32_t> anchor_positions = find_anchor(read_i_cigar,
                                                                        read_i_ref_offset,
                                                                        read_j_cigar,
                                                                        read_j_ref_offset);

            bool are_compatible = are_compatible_reads(read_i_seq, read_j_seq, anchor_positions, min_overlap_length);

            if (are_compatible) write_compatible_reads(read_i_id, read_j_id, output_basepath);

            // Limit the number of forward comparisons per read
            // This is very powerful in limiting search complexity for highly covered references
            // we go from quadratic O(n²) to amortised linear O(n*max_comparison_nb)
            if (++forward_comparison_count >= max_comparison_nb) break;
        }
    }
}

void clear_align_storage(single_ref_align_storage_t & align_storage)
{
    align_storage.ref_idx = -1;
    align_storage.read_ids.clear();
    align_storage.read_seqs.clear();
    align_storage.align_cigars.clear();
    align_storage.ref_offsets.clear();
    align_storage.read_offsets.clear();
}

void process_alignment_file(std::filesystem::path const & input_path,
                            reference_storage_t const & ref_storage,
                            std::filesystem::path const & output_basepath, 
                            std::uint16_t const min_overlap_length)
{
    debug_stream << "INFO: processing alignment file" << std::endl;
    time_t start_time = time(NULL);

    alignment_file_input input_file{input_path,
                                    ref_storage.ids,
                                    ref_storage.seqs,
                                    fields<field::ID,
                                           field::FLAG,
                                           field::REF_ID,
                                           field::REF_OFFSET,
                                           field::CIGAR,
                                           field::SEQ,
                                           field::OFFSET>{}};
    
    single_ref_align_storage_t align_storage{};

    for (auto & [read_id, flag, ref_idx, ref_offset, align_cigar, read_seq, read_offset] : input_file)
    {
        // debug_stream << "Stored ref_idx: " << align_storage.ref_idx << std::endl;

        if (not ref_idx) // Unmapped read
        {
            // print_debug(read_id, flag, ref_idx, ref_storage, ref_offset, align_cigar, read_seq, read_offset);
            continue;
        }
        else // Mapped read
        {
            if (ref_idx.value() != align_storage.ref_idx)
            {
                if (align_storage.ref_idx >= 0)
                {
                    find_reads_overlaps(align_storage, output_basepath, min_overlap_length);
                    clear_align_storage(align_storage);
                }
                align_storage.ref_idx = ref_idx.value();
            }

            // print_debug(read_id, flag, ref_idx, ref_storage, ref_offset, align_cigar, read_seq, read_offset);

            // If read was reverse-complemented before mapping, then revcomp its sequence
            if ((flag & 0x10) == 0x10) read_seq = read_seq | views::complement | std::views::reverse | views::to<std::vector>;
            
            // debug_stream << "final seq:   " << read_seq << std::endl;

            align_storage.read_ids.push_back(std::move(read_id));
            align_storage.read_seqs.push_back(std::move(read_seq)); // read_seq is already oriented as per the alignment strand
            align_storage.flags.push_back(std::move(flag));
            align_storage.align_cigars.push_back(std::move(align_cigar));
            align_storage.ref_offsets.push_back(ref_offset.value());
            align_storage.read_offsets.push_back(std::move(read_offset));
        }
    }

    if (align_storage.ref_idx >= 0 and align_storage.read_ids.size() > 1)
    {
        find_reads_overlaps(align_storage, output_basepath, min_overlap_length);
    }

    debug_stream << "INFO: finished processing alignment file " 
                 << "in " << time(NULL)-start_time << " seconds" << std::endl;
}

void run_program(std::filesystem::path const & input_path, 
                 std::filesystem::path const & reference_path,
                 std::filesystem::path const & output_basepath, 
                 std::uint16_t const min_overlap_length)
{
    // Read the reference file
    reference_storage_t ref_storage{};
    read_reference(reference_path, ref_storage);

    // Process the SAM/BAM file
    process_alignment_file(input_path, ref_storage, output_basepath, min_overlap_length);
}

struct cmd_arguments
{
    std::filesystem::path input_path{};
    std::filesystem::path reference_path{};
    std::filesystem::path output_basepath{"overlap_out"};
    uint16_t min_overlap_length{50};
};

void initialize_argument_parser(argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "Pierre Pericard";
    parser.info.short_description = "Compute MATAM overlap-graph from a SAM/BAM file";
    parser.info.version = "2.0.0b";
    parser.add_option(args.input_path, 'i', "input-sam", "Input SAM/BAM file",
                      option_spec::REQUIRED, input_file_validator{{"sam","bam"}});
    parser.add_option(args.reference_path, 'r', "reference", "Reference fasta file",
                      option_spec::REQUIRED, input_file_validator{{"fa","fasta"}});
    parser.add_option(args.output_basepath, 'o', "output", "Output basepath",
                      option_spec::DEFAULT);
    parser.add_option(args.min_overlap_length, 'm', "min-overlap", "Minimum overlap size",
                      option_spec::DEFAULT, arithmetic_range_validator{0, 65535});
}

int main(int argc, char ** argv)
{
    // argument_parser parser{"OvGraphBuild", argc, argv}; // initialise parser
    argument_parser parser{"OvGraphBuild", argc, argv, false}; // initialise parser // turn off version check
    // REMEMBER TO REACTIVATE VERSION CHECK
    cmd_arguments args{};
    
    initialize_argument_parser(parser, args);
    
    try
    {
         parser.parse();
    }
    catch (parser_invalid_argument const & ext)
    {
        std::cerr << "[Argument parsing error] " << ext.what() << "\n";
        return -1;
    }
    
    run_program(args.input_path, args.reference_path, args.output_basepath, args.min_overlap_length);
    
    return 0;
}