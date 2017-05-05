library(Rcpp)
delayedAssign(      "vcfGTz_num_SNPs", Rcpp:: cppFunction( plugins='cpp14', includes=c("#include<fstream>" ,'#include "/home/aaron/Research/Code/Imputation/ssimp_software/src/vcfGTz_reader.hh"')
	, 'IntegerVector vcfGTz_num_SNPs ( CharacterVector filename) {
			std:: string filename_ = {filename[0]};
			vcfGTz:: vcfGTz_reader reader(filename_);

			// read in the "magic" string
			auto magic = reader.read_string0();
			if(magic != "vcfGTz.0.0.2")
				stop("wrong file type?");

			// skip over the first block
			auto offset_over_this_block = reader.read_offset_at_start_of_block();
			reader.seek_relative(offset_over_this_block-8);

			// now read the start of the second block
			auto another_offset = reader.read_offset_at_start_of_block();
			auto description_of_second_block = reader.read_string0();
			if(description_of_second_block != "offsets.into.previous.block")
				stop("wrong file type?");

			int num_lines = reader.read_uint64_t();
			int num_SNPs = num_lines - 1; // as the first one is the header
			return {num_SNPs};
}'))
delayedAssign(      "vcfGTz_get_internal_offsets", Rcpp:: cppFunction( plugins='cpp14', includes=c("#include<fstream>" ,'#include "/home/aaron/Research/Code/Imputation/ssimp_software/src/vcfGTz_reader.hh"')
	, 'NumericVector vcfGTz_get_internal_offsets (CharacterVector filename, IntegerVector indices) {
			std:: string filename_ = {filename[0]};
			vcfGTz:: vcfGTz_reader reader(filename_);

			// read in the "magic" string
			auto magic = reader.read_string0();
			if(magic != "vcfGTz.0.0.2")
				stop("wrong file type?");

			// skip over the first block
			auto offset_over_this_block = reader.read_offset_at_start_of_block();
			reader.seek_relative(offset_over_this_block-8);

			// now read the start of the second block
			auto another_offset = reader.read_offset_at_start_of_block();
			auto description_of_second_block = reader.read_string0();
			if(description_of_second_block != "offsets.into.previous.block")
				stop("wrong file type?");

			int num_lines = reader.read_uint64_t();
			int num_SNPs = num_lines - 1; // as the first one is the header

			auto remember_this_position = reader.m_f.tellg();

			NumericVector results; // Numeric is like a 52-bit integer, in case the vcfGTz files are really big

			for( auto && i : indices ) {
				if(i>=1 && i<=num_SNPs) {
					reader.m_f.seekg( remember_this_position .operator+( 8*i ));
					// 0 would be the header, but we will use base-1
					// as this is R
					auto u64 = reader.read_uint64_t();
					results.push_back(u64);

					auto verify = results.at(results.size()-1);
					u64 == verify || (stop("problem with offsets? Really big file?"),false);
				}
				else
					results.push_back(NA_REAL);
			}

			return results;
}'))
delayedAssign(        "vcfGTz_get_1field_from_internal_offsets", Rcpp:: cppFunction( plugins='cpp14', includes=c("#include<fstream>" ,'#include "/home/aaron/Research/Code/Imputation/ssimp_software/src/vcfGTz_reader.hh"')
	, 'CharacterVector vcfGTz_get_1field_from_internal_offsets (CharacterVector filename, NumericVector file_offsets, IntegerVector which_field) {
			which_field.size() == 1 || (stop("only one field maybe be requested, between 1 and 7"),false);
			which_field.at(0) >= 1  || (stop("only one field maybe be requested, between 1 and 7"),false);
			which_field.at(0) <= 7  || (stop("only one field maybe be requested, between 1 and 7"),false);
			int number_of_columns_to_skip = which_field.at(0)-1;

			std:: string filename_ = {filename[0]};
			vcfGTz:: vcfGTz_reader reader(filename_);
			reader.m_f.good() || (stop(std::string("cannot find file: [") + filename_ + "]"),false);

			// read in the "magic" string
			auto magic = reader.read_string0();
			if(magic != "vcfGTz.0.0.2")
				stop("wrong file type?");

			// this time, we like the first block
			auto remember_this_position = reader.m_f.tellg();

			reader.read_offset_at_start_of_block(); // read this, but ignore it

			auto description_of_first_block = reader.read_string0();
			if(description_of_first_block != "manylines:GTonly:zlib")
				stop("wrong file type?");

			CharacterVector results;

			for(auto && o : file_offsets) {
				if(std::isnan(o)) {
					results.push_back(NA_STRING);
				}
				else {
					reader.m_f.seekg( remember_this_position.operator+( o ) );
					for(int i=0; i<number_of_columns_to_skip; ++i) {
						reader.read_smart_string0(); // to skip over fields
					}
					auto CHROM = reader.read_smart_string0();
					results.push_back(CHROM);
				}
			}

			return results;
}'))
delayedAssign(        "vcfGTz_get_012calls_from_internal_offsets", Rcpp:: cppFunction( plugins='cpp14', includes=c(
			    '#include "/home/aaron/Research/Code/Imputation/ssimp_software/src/R/ASSERT_for_R.hh"'
			,   "#include<fstream>"
			,   '#include "/home/aaron/Research/Code/Imputation/ssimp_software/src/vcfGTz_reader.hh"'
			)
	, 'IntegerMatrix vcfGTz_get_012calls_from_internal_offsets (CharacterVector filename, NumericVector file_offsets) {

			std:: string filename_ = {filename[0]};
			vcfGTz:: vcfGTz_reader reader(filename_);
			reader.m_f.good() || (stop(std::string("cannot find file: [") + filename_ + "]"),false);

			// read in the "magic" string
			auto magic = reader.read_string0();
			if(magic != "vcfGTz.0.0.2")
				stop("wrong file type?");

			// this time, we like the first block
			auto remember_this_position = reader.m_f.tellg();

			reader.read_offset_at_start_of_block(); // read this, but ignore it

			auto description_of_first_block = reader.read_string0();
			if(description_of_first_block != "manylines:GTonly:zlib")
				stop("wrong file type?");

			std:: cout << reader.read_smart_string0() << "\\n";
			std:: cout << reader.read_smart_string0() << "\\n";
			std:: cout << reader.read_smart_string0() << "\\n";
			std:: cout << reader.read_smart_string0() << "\\n";
			std:: cout << reader.read_smart_string0() << "\\n";
			std:: cout << reader.read_smart_string0() << "\\n";
			std:: cout << reader.read_smart_string0() << "\\n";

			auto doubly_compressed = reader.read_vector_of_char_with_leading_size();
			std:: cout << doubly_compressed.size() << "\\n";
			//auto un_z = zlib_vector:: inflate(doubly_compressed);

			IntegerMatrix results;
			results.push_back(0);
			results.push_back(1);
			results.push_back(2);
			results.push_back(3);
			//results.attr("dim") = Dimension(1,4);

			return results;
}'))
