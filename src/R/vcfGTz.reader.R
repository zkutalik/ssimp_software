library(Rcpp)
delayedAssign("vcfGTz_num_SNPs", Rcpp:: cppFunction( plugins='cpp14', includes=c("#include<fstream>" ,'#include "/home/aaron/Research/Code/Imputation/ssimp_software/src/vcfGTz_reader.hh"')
	, 'IntegerVector vcfGTz_num_SNPs( CharacterVector filename) {
			std:: string filename_ = {filename[0]};
			vcfGTz:: vcfGTz_reader reader(filename_);

			// read in the "magic" string
			auto magic = reader.read_string0();
			if(magic != "vcfGTz.0.0.2")
				return {-1};

			// skip over the first block
			auto offset_over_this_block = reader.read_offset_at_start_of_block();
			reader.seek_relative(offset_over_this_block-8);

			// now read the start of the second block
			auto another_offset = reader.read_offset_at_start_of_block();
			auto description_of_second_block = reader.read_string0();
			if(description_of_second_block != "offsets.into.previous.block")
				return {-1};

			int num_lines = reader.read_uint64_t();
			int num_SNPs = num_lines - 1; // as the first one is the header
			return {num_SNPs};
}'))
