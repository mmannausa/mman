#include <string>
using std::string;
#include <sstream>
using std::istringstream;
#include <vector>
using std::vector;
#include <unordered_map>
#include <algorithm>
using std::getline;
#include <iostream>
using std::endl;
using std::cout;
/*
	Matthew Mannausa
	Project 2
*/


bool IsValidDNASequence(const std::string& input) {
	char letter;
	for (long unsigned int i = 0; i < input.length(); ++i) {
		letter = input[i];//checks each letter of the sequence to see if it is valid
		if (letter != 'A' && letter != 'T' && letter != 'C' && letter != 'G') {
			return false;//it it is not valid returns false
		}
		
	}
	return true;// if each char passes then return true
}

string reverse_string(const string& input) {// a function just to return a reversed string
	string reverse;
	for (int i = input.size(); i >= 0; --i) {
		reverse += input[i];
	}
	return reverse;

}
void GetReverseComplementSequence(const std::string& input, std::string* const output) {
	string reversed_string;
	char letter;
	reversed_string = reverse_string(input);//reverses the string
	for (long unsigned int i = 0; i < reversed_string.size(); ++i) {// creates the anti-parrallel
		letter = reversed_string[i];
		if (letter == 'A') {
			*output += 'T';
		}
		else if (letter == 'T') {
			*output += 'A';
		}
		else if (letter == 'C') {
			*output += 'G';
		}
		else if (letter == 'G') {
			*output += 'C';
		}
	}

}
string rna_transcript(const string& input) {
	char letter;
	string RNA_string;
	for (long unsigned int i = 0; i < input.size(); ++i) {
		letter = input[i];
		if (letter == 'T') {//changes T's to U's
			RNA_string += 'U';
		}
		else {
			RNA_string += input[i];//if the letter is not a T just add it to the string
		}
	}
	return RNA_string;
}

std::string GetRNATranscript(const std::string& input) {
	string DNA_strand;
	string new_string, RNA_string;
	string reversed_string = reverse_string(input);
	char letter; 
	for (long unsigned int i = 0; i < reversed_string.size(); ++i) {// creates the anti-parrallel
		letter = reversed_string[i];
		if (letter == 'A') {
			DNA_strand += 'T';
		}
		else if (letter == 'T') {
			DNA_strand += 'A';
		}
		else if (letter == 'C') {
			DNA_strand += 'G';
		}
		else if (letter == 'G') {
			DNA_strand += 'C';
		}
	}
	
	RNA_string = rna_transcript(DNA_strand);
	return RNA_string;

}

std::vector<std::vector<std::string>> GetReadingFramesAsCodons(const std::string& input) {
	int x = 0;
	int y = 1;
	string rna_string = rna_transcript(input);
	string rna_tran = GetRNATranscript(input);
	vector<vector<string>> rna_vec(6); 
	string code; 

	for (int j = 0; j < 3; ++j) {
		x = 0;
		
		for (long unsigned int i = 0; rna_tran.size()/3 > i ; ++i) {
			
			code = rna_tran.substr(x, 3);//create the codons 
			rna_vec[j].push_back(code);
			x += 3;
		}
		
		rna_tran = rna_tran.substr(y);
		
	}
	
	for (int i = 3; i < 6; ++i) {
		x = 0;
		for (long unsigned int j = 0; j < rna_string.size()/3; ++j) {
			
			rna_vec[i].push_back(rna_string.substr(x, 3));//creates the codons
			x += 3;
		}
		
		rna_string = rna_string.substr(y);
	}
	return rna_vec;


	
}

std::string Translate(const std::vector<std::string>& codon_sequence) {
	char letter;
	string abriviation, vec_code;
	std::unordered_map<string, char> acid_map = {//unordered map from codes to corrisponding amino acids
		{"GCU",'A'},{"GCC",'A'},{"GCA",'A'},{"GCG",'A'},{"CGU",'R'},{"CGC",'R'},{"CGA",'R'},{"CGG",'R'},{"AGA",'R'},{"AGG",'R'},
		{"AAU",'N'},{"AAC",'N'},{"GAU",'D'},{"GAC",'D'},{"UGU",'C'},{"UGC",'C'},{"CAA",'Q'},{"CAG",'Q'},
		{"GAA",'E'},{"GAG",'E'},{"GGU",'G'},{"GGC",'G'},{"GGA",'G'},{"GGG",'G'},{"CAU",'H'},{"CAC",'H'},
		{"AUU",'I'},{"AUC",'I'},{"AUA",'I'},{"UUA",'L'},{"UUG",'L'},{"CUU",'L'},{"CUC",'L'},
		{"CUA",'L'},{"CUG",'L'},{"AAA",'K'},{"AAG",'K'},{"AUG",'M'},{"UUU",'F'},{"UUC",'F'},{"CCU",'P'},
		{"CCC",'P'},{"CCA",'P'},{"CCG",'P'},{"UCU",'S'},{"UCC",'S'},{"UCA",'S'},{"UCG",'S'},{"AGU",'S'},
		{"AGC",'S'},{"ACU",'T'},{"ACC",'T'},{"ACA",'T'},{"ACG",'T'},{"UGG",'W'},{"UAU",'Y'},{"UAC",'Y'},
		{"GUU",'V'},{"GUC",'V'},{"GUA",'V'},{"GUG",'V'},{"UAG",'*'},{"UGA",'*'},{"UAA",'*'} };
	

		for (long unsigned int j = 0; j < codon_sequence.size(); ++j) {
			vec_code = codon_sequence[j];
			letter = acid_map[vec_code];
			abriviation += letter;
			

		}

	
	return abriviation;




}

std::string GetLongestOpenReadingFrame(const std::string& DNA_sequence) {
	vector<vector<string>> dna_vec = GetReadingFramesAsCodons(DNA_sequence);
	vector<string> a_codon;
	string amino_acids;
	string long_read_frame, longest_read_frame = "", line;
	char first_letter;
	vector<string> read_frame_vec;
	
	for (long unsigned int j = 0; j < dna_vec.size(); ++j) {
		a_codon = dna_vec[j];//grab a single vector 
		amino_acids = Translate(a_codon); 
		istringstream read_frame(amino_acids);
		while(read_frame >> first_letter) {
			if (first_letter == 'M') {//reads until the first letter is a M
				long_read_frame = first_letter;
				while (read_frame >> first_letter) {
					long_read_frame += first_letter;//add the letters together
					if (first_letter == '*') {//stops if the letter is a *
						read_frame_vec.push_back(long_read_frame);

						if (long_read_frame.size() > longest_read_frame.size()) {//if the frame is longer then change the longest frame
							longest_read_frame = long_read_frame;

						}
						break;//break out of the inter while loop
					}

				}
			}
		}
	}
	
	return longest_read_frame;
}

