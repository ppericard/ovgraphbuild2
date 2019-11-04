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
    debug_stream << "read offset: " << read_offset << std::endl;
    debug_stream << "#################" << std::endl;
}

void find_reads_overlaps(single_ref_align_storage_t const & align_storage,
                         std::filesystem::path const & output_basepath, 
                         uint16_t const min_overlap_length)
{
    debug_stream << std::endl << std::endl;
    debug_stream << "Ref idx:  " << align_storage.ref_idx << std::endl;
    debug_stream << "Read names:  " << std::endl;
    debug_stream << align_storage.read_ids << std::endl;
    debug_stream << "#################" << std::endl << std::endl;

    size_t reads_number = align_storage.read_ids.size();

    if(reads_number>1)
    for(size_t i = 0; i < reads_number-1; ++i)
    {
        auto & read_i_id = align_storage.read_ids[i];
        auto & read_i_seq = align_storage.read_seqs[i];
        auto & read_i_cigar = align_storage.align_cigars[i];
        auto & read_i_ref_offset = align_storage.ref_offsets[i];

        for(size_t j = i+1; j < reads_number; ++j)
        {
            auto & read_j_id = align_storage.read_ids[j];
            auto & read_j_seq = align_storage.read_seqs[j];
            auto & read_j_cigar = align_storage.align_cigars[j];
            auto & read_j_ref_offset = align_storage.ref_offsets[j];

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
                            uint16_t const min_overlap_length)
{
    alignment_file_input input_file{input_path,
                                    ref_storage.ids,
                                    ref_storage.seqs,
                                    fields<field::ID,
                                           field::FLAG,
                                           field::REF_ID,
                                           field::REF_OFFSET,
                                           field::CIGAR,
                                           field::ALIGNMENT,
                                           field::SEQ,
                                           field::OFFSET>{}};
    
    single_ref_align_storage_t align_storage{};

    for (auto & [read_id, flag, ref_idx, ref_offset, align_cigar, alignment, read_seq, read_offset] : input_file)
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
            debug_stream << alignment << std::endl;

            // If read was reverse-complemented before mapping, then revcomp its sequence
            if ((flag & 0x10) == 0x10) read_seq = read_seq | views::complement | std::views::reverse | views::to<std::vector>;
            
            // debug_stream << "final seq:   " << read_seq << std::endl;

            align_storage.read_ids.push_back(std::move(read_id));
            align_storage.read_seqs.push_back(std::move(read_seq));
            align_storage.flags.push_back(std::move(flag));
            align_storage.align_cigars.push_back(std::move(align_cigar));
            align_storage.ref_offsets.push_back(ref_offset.value());
            align_storage.read_offsets.push_back(std::move(read_offset));
        }
    }

    if (align_storage.ref_idx >= 0)
    {
        find_reads_overlaps(align_storage, output_basepath, min_overlap_length);
    }
}

void run_program(std::filesystem::path const & input_path, 
                 std::filesystem::path const & reference_path,
                 std::filesystem::path const & output_basepath, 
                 uint16_t const min_overlap_length)
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
