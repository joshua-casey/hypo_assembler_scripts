#include <zlib.h> // for reading compressed .fq file
#include <iostream>
#include <ctime>
#include <chrono>
#include <unordered_set>
#include <unordered_map>
#include <tuple>
#include <vector>
#include <algorithm>
#include <htslib/kseq.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <omp.h>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/util.hpp>
using namespace std;
KSEQ_INIT(gzFile, gzread)

unsigned char seq_nt4_table[256] = {
        0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

int main(int argc, char* argv[]) {
    if(argc < 4) {
        cerr << "Usage: " << argv[0] << " <contigs fasta> <long reads alignments> <outputs>" << endl;
        return 1;
    
    }
    cerr << "Command: ";
    for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
    cerr << endl;
    
    vector<string> contig_names;
    vector<string> contigs;
    vector< vector<int> > supports;
    map<string, int> contig_name_to_id;
    
    cerr << "Reading contigs from " << argv[1] << endl;
    
    gzFile fp = gzopen(argv[1], "r");
    kseq_t *seq;
    seq = kseq_init(fp);
    int l;
    int ids = 0;
    while((l = kseq_read(seq)) >= 0) {
        std::string get_sequence = "";
        std::string get_id(seq->name.s, seq->name.l);
        
        for(int i = 0; i < seq->seq.l; i++) {
            get_sequence += seq->seq.s[i];
        }
        
        contigs.push_back(get_sequence);
        contig_names.push_back(get_id);
        
        supports.push_back(vector<int>(get_sequence.length()));
        
        contig_name_to_id[get_id] = ids;
        ids++;
    }
    
    
    cerr << "Reading bam file from " << argv[2] << endl;
    
    auto sam_file = sam_open(argv[2], "r");
    auto sam_header = sam_hdr_read(sam_file);
    auto current_align = bam_init1();
    
    while(sam_read1(sam_file, sam_header, current_align)>=0) {
        // ignore unmapped reads
        if (current_align->core.flag & (BAM_FUNMAP)) {
            continue;
        }
        
        // ignore secondary and failed reads
        if (current_align->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) { 
            continue;
        }
        
        // filter out quality value (long reads)
        if (current_align->core.qual < 2) {
            continue;
        }
        
        // filter out invalid reads
        if (current_align->core.pos < 0) {
            continue;
        }
        
        supports[current_align->core.tid][current_align->core.pos] += 1;
    }
    
    
    cerr << "Breaking misjoins." << endl;
    
    std::vector<std::string> new_contigs;
    std::vector<std::string> new_contig_ids;
    
    for(auto i = 0; i < contigs.size(); i++) {
        int a = contigs[i].size() * 3 / 4;
        int b = contigs[i].size() * 4 / 5;
        new_contig_ids.push_back(contig_names[i] + "_1");
        new_contig_ids.push_back(contig_names[i] + "_2");
        new_contigs.push_back(contigs[i].substr(0, b));
        new_contigs.push_back(contigs[i].substr(a));
    }
    
    
    cerr << "Writing misjoins to " << argv[3] << endl;
    
    std::ofstream ofile(argv[3]);
    if (!ofile.is_open()) {
        fprintf(stderr, "Error: File open error: Output File (%s) could not be opened!\n",argv[3]);
        exit(1);
    }
    for(int i = 0; i < new_contigs.size(); i++) {
        ofile << ">" << new_contig_ids[i] << "\n";
        ofile << new_contigs[i] << "\n";
    }
    ofile.close();
    
    return 0;
    
}
