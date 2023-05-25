#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <htslib/sam.h>
#include <vector>
#include "spoa/spoa.hpp"
#include "ssw.h"
#include "ssw_cpp.h"
#include <chrono>
#include <zlib.h>
extern "C" {
#include "kseq.h"
}
KSEQ_INIT(gzFile, gzread)
#include <unordered_map>
using std::string;
using std::cout;
using std::endl;
using std::istringstream;
using std::vector;
using std::getline;
using std::cerr;
using std::ifstream;
using std::unordered_map;
int usage();

unordered_map<string, string> load_reference_genome(const char *reference_file) {
    unordered_map<string, string> reference_genome;
    gzFile fp = gzopen(reference_file, "r");
    kseq_t *seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        reference_genome[seq->name.s] = seq->seq.s;
    }
    kseq_destroy(seq);
    gzclose(fp);
    return reference_genome;
}

//g++ STR.cpp -o STR -I../../htslib/include -L../../htslib/lib -lhts -I../include -L../lib -lspoa -lz -Wall ssw_cpp.cpp ssw.c

static void PrintAlignment(const StripedSmithWaterman::Alignment& alignment, const string& ref, const string& query, const int& length , const int& coverage, const string& chrom, const int& start, const int& end, const int& flanking_len,bool isOutputOther ) {
	if (isOutputOther) {
		cout << endl;
		cout << "===== SSW result =====" << endl;
		cout << "Best Smith-Waterman score:\t" << alignment.sw_score << endl
	     	     << "Next-best Smith-Waterman score:\t" << alignment.sw_score_next_best << endl
	             << "Reference start:\t" << alignment.ref_begin << endl
	             << "Reference end:\t" << alignment.ref_end << endl
	             << "Query start:\t" << alignment.query_begin << endl
	             << "Query end:\t" << alignment.query_end << endl
	             << "Next-best reference end:\t" << alignment.ref_end_next_best << endl
	             << "Number of mismatches:\t" << alignment.mismatches << endl
	             << "Cigar: " << alignment.cigar_string << endl;
	}

	int ref_pos = alignment.ref_begin;
	//int query_pos = alignment.query_begin;
	int query_pos = 0;
	int lengthvariation=0;
	int site=flanking_len;
	site--;

	for (const auto& c : alignment.cigar) {
		char op = c & 0xf;
		uint32_t len = c >> 4;

		if (op == 1) {  // Insertion in the query
			if (isOutputOther) {
				cout << "Insertion in query: " << query.substr(query_pos, len) << " at reference position " << ref_pos << endl;
			}
			if( ref_pos >= site && ref_pos <= site+length ){
				lengthvariation+=len;
			}
			query_pos += len;
		}
		else if (op == 2) {  // Deletion in the query
			if (isOutputOther) {
				cout << "Deletion in query: " << ref.substr(ref_pos, len) << " at reference position " << ref_pos << endl;
			}
			if( ref_pos >= site && ref_pos <= site+length ){
				lengthvariation-=len;
			}
			ref_pos += len;
		}
		else if (op == 4) {  // Soft clipping in the query
			query_pos += len;
		}
		else if (op == 7 || op == 0 || op == 8 ) {  // Match in both ('=' or 'M')
			ref_pos += len;
			query_pos += len;
		}
	}

	if( isOutputOther ){
		if( coverage >=3 ){
			cout << "Result: " << chrom << "\t" << start << "\t" << end << "\t" << coverage << "\t" << lengthvariation << endl;
		}
		else{
			cout << "Result: " << chrom << "\t" << start << "\t" << end << "\t" << coverage << "\t" << "NA" << endl;
		}
	}
	else{
		if( coverage >=3 ){
                        cout << chrom << "\t" << start << "\t" << end << "\t" << coverage << "\t" << lengthvariation << endl;
                }
                else{
                        cout << chrom << "\t" << start << "\t" << end << "\t" << coverage << "\t" << "NA" << endl;
                }
	}

	if (isOutputOther) {
		cout << "======================" << endl;
	}
}


void process_bed_line(const string &line, samFile *in, bam_hdr_t *header, const char *alignment_file, int min_mapping_quality, unordered_map<string, string>& reference_genome, bool isOutputOther) {

	istringstream iss(line);
	string chrom;
	int start;
	int end;
	iss >> chrom >> start >> end;
	int str_len=end-start+1;
	//cout << line << endl;

	// load index file
	hts_idx_t *idx = sam_index_load(in, alignment_file);
	if (idx == NULL) {
		cerr << "Error loading index for BAM file" << endl;
		return;
	}

	int tid = bam_name2id(header, chrom.c_str());
	hts_itr_t *iter = sam_itr_queryi(idx, tid, start, end);

	bam1_t *aln = bam_init1();
	vector<string> all_reads;//a vector to stroge all the qualified reads
	int coverage=0;

	if (isOutputOther) {
		cout << line << endl;
		cout << "===reads sequences in bam file===" << endl;
	}
	while (sam_itr_next(in, iter, aln) >= 0) { //process each reads

		// Filter out alignments with mapping quality < min_mapping_quality
		if (aln->core.qual < min_mapping_quality) {
			continue;
		}

		uint32_t *cigar = bam_get_cigar(aln);
		int n_cigar = aln->core.n_cigar;
		bool has_clipping_at_ends = false;

		int clipping_cutoff;
		uint8_t *seq = bam_get_seq(aln);
		int seq_length = aln->core.l_qseq;
		if( seq_length <= 100 ){
			clipping_cutoff=35;
		}
		else if( seq_length <= 150 ){
			clipping_cutoff=50;
		}
		else if( seq_length <= 200 ){
			clipping_cutoff=75;
		}

		// Check for soft or hard clipping at the beginning
		int first_op = bam_cigar_op(cigar[0]);
		int first_op_len = bam_cigar_oplen(cigar[0]);
		if ((first_op == BAM_CSOFT_CLIP || first_op == BAM_CHARD_CLIP) && first_op_len > clipping_cutoff ) {
			has_clipping_at_ends = true;
		}

		// Check for soft or hard clipping at the end
		int last_op = bam_cigar_op(cigar[n_cigar - 1]);
		int last_op_len = bam_cigar_oplen(cigar[n_cigar - 1]);
		if ((last_op == BAM_CSOFT_CLIP || last_op == BAM_CHARD_CLIP) && last_op_len > clipping_cutoff ) {
			has_clipping_at_ends = true;
		}

		// Skip the reads if it has soft or hard clipping at both ends with lenght > 50;
		if (has_clipping_at_ends) {
			continue;
		}

		//uint8_t *seq = bam_get_seq(aln);
		//int seq_len = aln->core.l_qseq;
		string reads;
		coverage ++;

		// print the qualified reads
		if (isOutputOther) {
			cout << bam_get_qname(aln) << ": ";
		}
		for (int i = 0; i < seq_length; ++i) {
			char base = seq_nt16_str[bam_seqi(seq, i)];
			reads += base;
			if (isOutputOther) {
				cout << base;
			}
		}
		if (isOutputOther) {
			cout << endl;
		}
		all_reads.push_back(reads);
	}

	if(isOutputOther){

		if( coverage == 0 ){
			cout << "Result: " << chrom << "\t" << start << "\t" << end << "\t" << coverage << "\t" << "NA" << endl;
			return;
		}
		else if( coverage > 200 ){
			cout << "Result: " << chrom << "\t" << start << "\t" << end << "\t" << coverage << "\t" << "NA" << endl;
			return;
		}
	}
	else{
		if( coverage == 0 ){
			cout << chrom << "\t" << start << "\t" << end << "\t" << coverage << "\t" << "NA" << endl;
			return;
		}
		else if( coverage > 200 ){
			cout << chrom << "\t" << start << "\t" << end << "\t" << coverage << "\t" << "NA" << endl;
			return;
		}
	}

	//Next, spoa will be used to generate the consensus sequence of the read
	//Pay attention, currently, only one consensus sequence will be created. So it only works well on highly homozygous genome, like Arabidopsis thaliana.
	
	//kNW for Needleman-Wunsch (global alignment)
	//kSW for Smith-Waterman (local alignment)
	//kOV for overlap (semi-global alignment)

	auto alignment_engine = spoa::AlignmentEngine::Create(
	spoa::AlignmentType::kOV, 5, -4, -8, -6, -8, -6);

	spoa::Graph graph{};
	for (const auto& it : all_reads) {
		auto alignment = alignment_engine->Align(it, graph);
		graph.AddAlignment(alignment, it);
	}

	auto consensus = graph.GenerateConsensus();
	size_t consensus_length = consensus.length();//problem !!!
	//cout << "Length" << consensus_length << endl;
	
	if (isOutputOther) {
		cout << endl;
		cout << "===Multiple sequence alignment===" << endl;
	}
	auto msa = graph.GenerateMultipleSequenceAlignment();

	if (isOutputOther) {
		for (const auto& it : msa) {
			cout << it << endl;
		}
		cout << endl;

		cout << "===Coverage:" << "\t" << coverage << endl;
		cout << endl;
	}

	//To make a decent alignment, here I will extract a longer reference sequence
	int flanking_len;
	if( consensus_length <= 100 ){
		flanking_len=60;
	}
	else if( consensus_length <= 200 ){
	       flanking_len=110;
	}
	else if( consensus_length <= 300 ){
		flanking_len=160;
	}
	else if( consensus_length <= 400 ){
		flanking_len=210;
	}
	else if( consensus_length <= 500 ){
		flanking_len=260;
	}
	else{
		flanking_len=500;
	}

	//Then, consensus sequence will be aligned to the corresponding region on the reference genome using ssw.
	string reference_seq = reference_genome.at(chrom).substr(std::max(0, start - flanking_len), end - start + flanking_len*2);

	int match_score = 2;
	int mismatch_penalty = 4;
	int gap_open_penalty = 10;
	int gap_extend_penalty = 1;

        int32_t maskLen = strlen(consensus.c_str())/2;
        maskLen = maskLen < 50 ? 50 : maskLen;
        // Declares a default Aligner
        StripedSmithWaterman::Aligner aligner(match_score, mismatch_penalty, gap_open_penalty, gap_extend_penalty);
	//StripedSmithWaterman::Aligner aligner;
	// Declares a default filter
        StripedSmithWaterman::Filter filter;
        // Declares an alignment that stores the result
        StripedSmithWaterman::Alignment alignment;
        // Aligns the query to the ref
        aligner.Align(consensus.c_str(), reference_seq.c_str(), reference_seq.size(), filter, &alignment, maskLen);
	if (isOutputOther) {
		cout << "================" << endl;
       		cout << ">Consensus" << endl;
		cout << consensus << endl;
		cout << ">Reference_seq" << endl;
        	cout << reference_seq << endl;
		cout << "================" << endl;
		cout << endl;
	}
        PrintAlignment(alignment, reference_seq, consensus, str_len, coverage, chrom, start, end, flanking_len, isOutputOther);
	//cout << endl;

	bam_destroy1(aln);
	hts_itr_destroy(iter);
	hts_idx_destroy(idx);
}

int main(int argc, char **argv) {
	int parameters_num=0;
	int parameters;
	int min_mapping_quality=10;
	bool isOutputOther=false;
	const char *bed_file;
	const char *reference_file;
	const char *alignment_file;

	while(( parameters=getopt(argc, argv, "a:b:f:m:nh")) >= 0){
		if( parameters == 'a' ){
			alignment_file = optarg;
			parameters_num++;
		}
		else if( parameters == 'b' ){
			bed_file = optarg;
			parameters_num++;
		}
		else if( parameters == 'f' ){
			reference_file = optarg;
			parameters_num++;
		}
		else if( parameters == 'm' ){
			min_mapping_quality = atof(optarg);
			parameters_num++;
		}
		else if( parameters == 'n' ){
			isOutputOther = true;
			parameters_num++;
		}
		else if( parameters == 'h' ){
			return usage();
		}
	}

	if( parameters_num == 0 ){
		return usage();
	}
	else if( alignment_file == nullptr || bed_file == nullptr || reference_file == nullptr){
		cerr << "Your should input at least three parameters, -a BAM/CRAM -b BED -f Reference" << endl;
		return 1;
	}

	// open the alignment file
	samFile *sam_file = sam_open(alignment_file, "r");
	if (sam_file == nullptr) {
		cerr << "Error: Failed to open alignment file" << endl;
		return 1;
	}

	// Set the reference genome if it's a CRAM file
	if (hts_get_format(sam_file)->format == htsExactFormat::cram) {
		if (hts_set_fai_filename(sam_file, reference_file) != 0) {
			cerr << "Error: Failed to set reference genome for CRAM file" << endl;
			sam_close(sam_file);
			return 1;
		}
	}

	// check the alignment file
	bam_hdr_t *bam_header = sam_hdr_read(sam_file);
	if (bam_header == nullptr) {
		cerr << "Error: Failed to read header" << endl;
		sam_close(sam_file);
		return 1;
	}

	// open the bed file
	ifstream bed_stream(bed_file);
	if (!bed_stream) {
		cerr << "Error: Failed to open BED file" << endl;
		bam_hdr_destroy(bam_header);
		sam_close(sam_file);
		return 1;
	}

	// load genome to unorder_map Chr->sequence;
	auto reference_genome = load_reference_genome(reference_file);

	//loop through all the records in the bed file
	string line;
	while (getline(bed_stream, line)) {
		//record the time
		auto start_time = std::chrono::high_resolution_clock::now();
		//process each record
		process_bed_line(line, sam_file, bam_header, alignment_file, min_mapping_quality, reference_genome, isOutputOther);
		auto end_time = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
		if(isOutputOther){
			cout << "===Time taken for this iteration: " << duration.count() << " milliseconds===" << endl;
			cout << endl;
		}
	}

	bed_stream.close();
	bam_hdr_destroy(bam_header);
	sam_close(sam_file);

	return 0;
}

int usage(){
	cerr << "Usage:" << "genoTR BAM/CRAM BED Reference" << endl;
	cerr << "arguments:" << endl;
	cerr << "       -a string     alignment file" << endl;
	cerr << "       -b string     bed file, should not include header" << endl;
	cerr << "       -f string     reference file" << endl;
	cerr << "       -m int        minimum cutoff for mapping quality" << endl;
        cerr << "       -n            noisy but informatics output, that help you know how this program work" << endl;
	cerr << "       -h            help" << endl;
	return 1;
}
