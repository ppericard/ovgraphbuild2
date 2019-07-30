#include <seqan3/argument_parser/all.hpp>     // includes all necessary headers
#include <seqan3/core/debug_stream.hpp>  // our custom output stream
#include <seqan3/std/filesystem>          // use std::filesystem::path

using namespace seqan3;

struct cmd_arguments
{
    std::filesystem::path input_path{};
	std::filesystem::path output_basepath{};
    uint16_t min_overlap_length{50};
};

void initialize_argument_parser(argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "Pierre Pericard";
    parser.info.short_description = "Compute MATAM overlap-graph from a SAM/BAM file";
    parser.info.version = "2.0.0b";
	
	parser.add_option(args.input_path, 'i', "input", "Input SAM/BAM file",
					  option_spec::REQUIRED);
	parser.add_option(args.output_basepath, 'o', "output", "Output basepath",
					  option_spec::REQUIRED);
	parser.add_option(args.min_overlap_length, 'm', "min-overlap", "Minimum overlap size",
					  option_spec::DEFAULT, arithmetic_range_validator{0, 65535});
}

int main(int argc, char ** argv)
{
    argument_parser myparser{"OvGraphBuild", argc, argv};        // initialise myparser
    cmd_arguments args{};
    initialize_argument_parser(myparser, args);
    try
    {
         myparser.parse();                                          // trigger command line parsing
    }
    catch (parser_invalid_argument const & ext)                     // catch user errors
    {
        debug_stream << "[Argument parsing error] " << ext.what() << "\n"; // customise your error message
        return -1;
    }
    // parsing was successful !
    // we can start running our progra
	// run_ovgraphbuild();
	
    debug_stream << "Hello world\n";
    return 0;
}
