//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iostream>
#include <random>
#include <algorithm>
#include <zlib.h>

#include "sequence_container.h"
#include "logger.h"

#include "palisade_header.h"

using namespace std;


size_t SequenceContainer::g_nextSeqId = 0;

const FastaRecord::Id FastaRecord::ID_NONE = 
			Id(std::numeric_limits<uint32_t>::max());


bool SequenceContainer::isFasta(const std::string& fileName)
{
	std::string withoutGz = fileName;
	if (fileName.substr(fileName.size() - 3) == ".gz")
	{
		withoutGz = fileName.substr(0, fileName.size() - 3);
	}

	size_t dotPos = withoutGz.rfind(".");
	if (dotPos == std::string::npos)
	{
		throw ParseException("Can't identify input file type");
	}
	std::string suffix = withoutGz.substr(dotPos + 1);

	if (suffix == "fasta" || suffix == "fa")
	{
		return true;
	}
	else if (suffix == "fastq" || suffix == "fq")
	{
		return false;
	}
	throw ParseException("Can't identify input file type");
}

// Add Fasta record
// Fill std::unordered_map<std::string, FastaRecord::Id> _nameIndex and vector<FastaRecord> _seqIndex
FastaRecord::Id SequenceContainer::addSequence(const FastaRecord& seqRec)
{
	if (!_offsetInitialized)
	{
		_offsetInitialized = true;
		_seqIdOffest = g_nextSeqId;
	}
	FastaRecord::Id newId(g_nextSeqId);
	if (_seqIndex.size() != g_nextSeqId - _seqIdOffest) 
	{
		throw std::runtime_error("something wrong with sequence ids!");
	}
	g_nextSeqId += 2;

	cout<<"newId: "<<newId<<"; seqRec.description: "<<seqRec.description<<endl;

	_seqIndex.emplace_back(seqRec.sequence, "+" + seqRec.description, 		// SequenceIndex = vector<FastaRecord>
						   newId);										// Append current Fasta record to vector<FastaRecord> _seqIndex

	if (_nameIndex.count(_seqIndex.back().description))
	{
		throw ParseException("The input contain reads with duplicated IDs. "
							 "Make sure all reads have unique IDs and restart. "
							 "The first problematic ID was: " +
			 				 _seqIndex.back().description.substr(1));
	}

	// FastaRecord class has fields 'string description' and 'class Id id'; hash map _nameIndex used to map string to id
	_nameIndex[_seqIndex.back().description] = _seqIndex.back().id;

	_seqIndex.emplace_back(seqRec.sequence.complement(), 						// reverse complement-> A-T, G-C
						   "-" + seqRec.description, newId.rc());				// newId starts with 0, then, newId.rc()=1;
	_nameIndex[_seqIndex.back().description] = _seqIndex.back().id;

	return _seqIndex.back().id.rc();
}

void SequenceContainer::loadFromFile(const std::string& fileName, 
									 int minReadLength)
{
	std::vector<FastaRecord> records;
	if (this->isFasta(fileName))
	{
		this->readFasta(records, fileName);				// read all FASTA into vector records
	}
	else
	{
		this->readFastq(records, fileName);
	}

	cout<<"records.size(): "<<records.size()<<endl;
	
	//shuffling input reads
	//std::vector<size_t> indicesPerm(records.size());
	//for (size_t i = 0; i < indicesPerm.size(); ++i) indicesPerm[i] = i;
	//std::random_shuffle(indicesPerm.begin(), indicesPerm.end());

	// Append each FASTA recrod to vector<FastaRecord> _seqIndex via addSequence function
	for (size_t i = 0; i < records.size(); ++i)
	{
		if (records[i].sequence.length() > (size_t)minReadLength)
		{
			this->addSequence(records[i]);			// Each Fasta record with sequence length > minReadLength is added to SequenceContainer
		}
	}
}

int SequenceContainer::computeNxStat(float fraction) const
{
	std::vector<int32_t> readLengths;
	int64_t totalLengh = 0;
	for (const auto& read : _seqIndex) 
	{
		readLengths.push_back(read.sequence.length());
		totalLengh += read.sequence.length();
	}
	std::sort(readLengths.begin(), readLengths.end(),
			  [](int32_t a, int32_t b) {return a > b;});

	int32_t nx = 0;
    int64_t cummulativeLen = 0;
	for (auto l : readLengths)
	{
        cummulativeLen += l;
        if (cummulativeLen > fraction * totalLengh)
		{
            nx = l;
            break;
		}
	}
	return nx;
}

//adds sequence ad it's complement
const FastaRecord& 
	SequenceContainer::addSequence(const DnaSequence& sequence, 
								   const std::string& description)
{
	auto newId = this->addSequence({sequence, description, 
								   FastaRecord::ID_NONE});
	return _seqIndex[newId._id - _seqIdOffest];
}

// vector<FastaRecord> records;
// struct FastaRecord has fields - Id id; DnaSequence sequence;
size_t SequenceContainer::readFasta(std::vector<FastaRecord>& record, 
									const std::string& fileName)
{
	size_t BUF_SIZE = 32 * 1024 * 1024;
	char* rawBuffer = new char[BUF_SIZE];
	auto* fd = gzopen(fileName.c_str(), "rb");
	if (!fd)
	{
		throw ParseException("Can't open reads file");
	}

	record.clear();
	int lineNo = 1;
	std::string header; vecCT header_ct;
	std::string sequence; vecCT sequence_ct;
	std::string nextLine; vecCT nextLine_ct;

	try
	{
		while(!gzeof(fd))
		{
			// cout<<"################################################"<<endl;
			
			//get a new line
			for (;;)
			{
				char* read = gzgets(fd, rawBuffer, BUF_SIZE);			// read data from compressed fd till \n hit or BUF_SIZE-1 chars read
				if (!read) break;
				
				// nextLine += read;										// append read to line;  
				// std::cout<<"nextLine: "<<nextLine<<std::endl;
				// if (nextLine.empty()) break;

				int orig_size=nextLine_ct.size();
				// printf("orig_size: %d\n",orig_size);
				// printf("strlen(read): %d\n",(int)(strlen(read)));
				nextLine_ct.resize(orig_size+ceil(strlen(read)*1.0/16384));	vecInt v(16384,0);
				int i;
				for(i=0;i<strlen(read);i++){
					if(i!=0 && i%16384==0){
						nextLine_ct[orig_size+i/16384-1]=encrypt_plaintext_vector_to_ciphertext(v);
						std::fill(v.begin(),v.end(),0);
					}
					v[i%16384]=read[i];
				}
				// cout<<"End i: "<<i<<endl;
				nextLine_ct[orig_size+(i-1)/16384]=encrypt_plaintext_vector_to_ciphertext(v);
				// cout<<"nextLine_ct.size(): "<<nextLine_ct.size()<<endl;

				if(nextLine_ct.size()==0||(nextLine_ct.size()==1 && decrypt_ciphertext_to_plaintext_vector(nextLine_ct[0])[0]==0))	break;

				// if (nextLine.back() == '\n')							// If read ended with \n, pop out the \n and break- we have read one line
				// {
				// 	nextLine.pop_back();
				// 	cout<<"break1"<<endl;
				// 	// break;
				// }

				CT c1=nextLine_ct[orig_size+(i-1)/16384];
				vecInt c1_vec=decrypt_ciphertext_to_plaintext_vector(c1);
				
				if(c1_vec[(i-1)%16384]=='\n'){
					c1_vec[(i-1)%16384]=0;
					nextLine_ct[orig_size+(i-1)/16384]=encrypt_plaintext_vector_to_ciphertext(c1_vec);
					// cout<<"break2"<<endl;
					break;
				}
				
			}

			// if (nextLine.empty()) continue;
			if(nextLine_ct.size()==0||(nextLine_ct.size()==1 && decrypt_ciphertext_to_plaintext_vector(nextLine_ct[0])[0]==0))	continue;

			// if (nextLine.back() == '\r') nextLine.pop_back();
			
			vecInt d=decrypt_ciphertext_to_plaintext_vector(nextLine_ct[nextLine_ct.size()-1]);
			int bck=0;
			for(bck=0;bck<16384;bck++){
				if(d[bck]==0)
					break;
			}
			bck--;
			// printf("nextLine.back(): %c, d[bck]: %c\n",nextLine.back(),d[bck]);

			if(d[bck]=='\r')	d[bck]=0;
			nextLine_ct[nextLine_ct.size()-1]=encrypt_plaintext_vector_to_ciphertext(d);

			// /cout<<"nextLine[0]: "<<(char)(nextLine[0])<<"; decrypt_ciphertext_to_plaintext_vector(nextLine_ct[0])[0]: "<<(char)(decrypt_ciphertext_to_plaintext_vector(nextLine_ct[0])[0])<<endl;

			// if (nextLine[0] == '>')
			if(decrypt_ciphertext_to_plaintext_vector(nextLine_ct[0])[0]=='>')
			{
				// if (!header.empty())
				// out<<"nextLine is a header\n"<<endl;
				if(!(header_ct.size()==0||(header_ct.size()==1 && decrypt_ciphertext_to_plaintext_vector(header_ct[0])[0]==0)))
				{
					// cout<<"Header is not empty"<<endl;
					if (sequence_ct.size()==0||(sequence_ct.size()==1 && decrypt_ciphertext_to_plaintext_vector(sequence_ct[0])[0]==0)) throw ParseException("empty sequence");

					// cout<<"Adding sequence: "<<sequence<<"\n with header: "<<header<<endl;
					string sequence_ct_string=""; string header_ct_string="";
					for(int i=0;i<sequence_ct.size();i++){
						vecInt vv=decrypt_ciphertext_to_plaintext_vector(sequence_ct[i]);
						for(int j=0;j<vv.size();j++)
						{
							if(vv[j]==0)
								break;
							sequence_ct_string+=(char)(vv[j]);
						}
					}
					for(int i=0;i<header_ct.size();i++){
						vecInt vv=decrypt_ciphertext_to_plaintext_vector(header_ct[i]);
						for(int j=0;j<vv.size();j++)
						{
							if(vv[j]==0)
								break;
							header_ct_string+=(char)(vv[j]);
						}
					}
					
					// cout<<"sequence_ct_string: "<<sequence_ct_string<<"\n header_ct_string: "<<header_ct_string<<endl;

					record.emplace_back(DnaSequence(sequence_ct_string), header_ct_string, 				// nextLine currently has new header; old header and corresponding sequence are added to record
										FastaRecord::ID_NONE);
					// cout<<"record.back().id: "<<record.back().id<<"; record.back().description: "<<record.back().description<<endl;
					// for(int j=0;j<record.size();j++){
					// 	cout<<"record["<<j<<"].id: "<<record[j].id<<"; ";
					// }
					// cout<<endl;

					// sequence.clear();
					// header.clear();

					sequence_ct.clear(); header_ct.clear();
				}

				// this->validateHeader(nextLine);
				this->validateHeader_ct(nextLine_ct);

				// header = nextLine;
				header_ct=nextLine_ct;
			}
			else
			{
				// cout<<"validate sequence with nextLine"<<endl;
				// this->validateSequence(nextLine);
				this->validateSequence_ct(nextLine_ct);					// header for this sequence has already been read in previous iteration
				// std::copy(nextLine.begin(), nextLine.end(), 	// current sequence, present in nextLine, is copied to variable sequence
				// 		  std::back_inserter(sequence));
				std::copy(nextLine_ct.begin(), nextLine_ct.end(), 	
						  std::back_inserter(sequence_ct));
			}

			++lineNo;
			// nextLine.clear();
			nextLine_ct.clear();

		}
		
		// if (sequence.empty()) throw ParseException("empty sequence");
		if (sequence_ct.size()==0||(sequence_ct.size()==1 && decrypt_ciphertext_to_plaintext_vector(sequence_ct[0])[0]==0)) throw ParseException("empty sequence");

		// if (header.empty())
		if (header_ct.size()==0||(header_ct.size()==1 && decrypt_ciphertext_to_plaintext_vector(header_ct[0])[0]==0))
		{
			throw ParseException("Fasta fromat error");
		}

		string sequence_ct_string=""; string header_ct_string="";
		for(int i=0;i<sequence_ct.size();i++){
			vecInt vv=decrypt_ciphertext_to_plaintext_vector(sequence_ct[i]);
			for(int j=0;j<vv.size();j++)
			{
				if(vv[j]==0)
					break;
				sequence_ct_string+=(char)(vv[j]);
			}
		}
		for(int i=0;i<header_ct.size();i++){
			vecInt vv=decrypt_ciphertext_to_plaintext_vector(header_ct[i]);
			for(int j=0;j<vv.size();j++)
			{
				if(vv[j]==0)
					break;
				header_ct_string+=(char)(vv[j]);
			}
		}

		record.emplace_back(DnaSequence(sequence_ct_string), header_ct_string, 
							FastaRecord::ID_NONE);

	}
	
	catch (ParseException& e)
	{
		std::stringstream ss;
		ss << "parse error in " << fileName << " on line " << lineNo << ": " << e.what();
		gzclose(fd);
		throw ParseException(ss.str());
	}

	delete[] rawBuffer;
	gzclose(fd);
	return record.size();
}

size_t SequenceContainer::readFastq(std::vector<FastaRecord>& record, 
									const std::string& fileName)
{

	size_t BUF_SIZE = 32 * 1024 * 1024;
	char* rawBuffer = new char[BUF_SIZE];
	auto* fd = gzopen(fileName.c_str(), "rb");
	if (!fd)
	{
		throw ParseException("Can't open reads file");
	}

	record.clear();
	int lineNo = 1;
	int stateCounter = 0;
	std::string header; 
	std::string nextLine;
	try
	{
		while (!gzeof(fd))
		{
			//get a new line
			for (;;)
			{
				char* read = gzgets(fd, rawBuffer, BUF_SIZE);
				if (!read) break;
				nextLine += read;
				if (nextLine.empty()) break;
				if (nextLine.back() == '\n')
				{
					nextLine.pop_back();
					break;
				}
			}

			if (nextLine.empty()) 
			{
				stateCounter = (stateCounter + 1) % 4;
				continue;
			}
			if (nextLine.back() == '\r') nextLine.pop_back();

			if (stateCounter == 0)
			{
				if (nextLine[0] != '@') throw ParseException("Fastq format error");
				header = nextLine;
				this->validateHeader(header);
			}
			else if (stateCounter == 1)
			{
				this->validateSequence(nextLine);
				record.emplace_back(DnaSequence(nextLine), header, 
									FastaRecord::ID_NONE);
			}
			else if (stateCounter == 2)
			{
				if (nextLine[0] != '+') throw ParseException("Fastq fromat error");
			}
			stateCounter = (stateCounter + 1) % 4;
			++lineNo;
			nextLine.clear();
		}
	}
	catch (ParseException& e)
	{
		std::stringstream ss;
		ss << "parse error in " << fileName << " on line " << lineNo << ": " << e.what();
		gzclose(fd);
		throw ParseException(ss.str());
	}

	gzclose(fd);
	delete[] rawBuffer;
	return record.size();
}


void SequenceContainer::validateHeader(std::string& header)
{
	size_t delim = 0;
	for (delim = 0; delim < header.length(); ++delim)
	{
		if (std::isspace(header[delim])) break;
	}

	header = header.substr(1, delim - 1);
	if (header.empty()) throw ParseException("empty header");
}

void SequenceContainer::validateHeader_ct(vecCT& header_ct)
{
	int flag=0;

	for(int i=0;i<header_ct.size();i++){
		vecInt vv=decrypt_ciphertext_to_plaintext_vector(header_ct[i]);
		for(int j=0;j<vv.size();j++)
		{
			if(vv[j]==' ')
			{
				flag=1;
				for(int k=j;k<vv.size();k++)
					vv[k]=0;
				header_ct[i]=encrypt_plaintext_vector_to_ciphertext(vv);

				for(int k=i+1;k<header_ct.size();k++)
				{
					vecInt vv2=decrypt_ciphertext_to_plaintext_vector(header_ct[k]);
					fill(vv2.begin(),vv2.end(),0);
					header_ct[k]=encrypt_plaintext_vector_to_ciphertext(vv2);
				}

				break;
			}
		}

		if(flag==1)
			break;
	}

	vecInt vv=decrypt_ciphertext_to_plaintext_vector(header_ct[0]);
	vecInt vv2(16384,0);
	copy(vv.begin()+1,vv.end(),vv2.begin());
	header_ct[0]=encrypt_plaintext_vector_to_ciphertext(vv2);

	if (header_ct.size()==0 || (header_ct.size()==1 && decrypt_ciphertext_to_plaintext_vector(header_ct[0])[0]==0)) throw ParseException("empty header");
}

void SequenceContainer::validateSequence(std::string& sequence)
{
	const std::string VALID_CHARS = "ACGT";
	for (size_t i = 0; i < sequence.length(); ++i)
	{
		if (DnaSequence::dnaToId(sequence[i]) == -1U)
		{
			sequence[i] = VALID_CHARS[rand() % 4];
		}
	}
}

void SequenceContainer::validateSequence_ct(vecCT& sequence_ct)
{
	const std::string VALID_CHARS = "ACGT";

	for (size_t i = 0; i < sequence_ct.size(); ++i)
	{
		vecInt vv=decrypt_ciphertext_to_plaintext_vector(sequence_ct[i]); int flag=0;

		for(int j=0;j<vv.size();j++){
			if(vv[j]!=0){
				if (DnaSequence::dnaToId(vv[j]) == -1U)
				{
					flag=1;
					vv[j] = VALID_CHARS[rand() % 4];
				}
			}
		}

		if(flag==1)
			sequence_ct[i]=encrypt_plaintext_vector_to_ciphertext(vv);
	}
}

void SequenceContainer::writeFasta(const std::vector<FastaRecord>& records, 
								   const std::string& filename,
								   bool onlyPositiveStrand)
{
	static const size_t FASTA_SLICE = 80;

	Logger::get().debug() << "Writing FASTA";
	FILE* fout = fopen(filename.c_str(), "w");
	if (!fout) throw std::runtime_error("Can't open " + filename);
	
	for (const auto& rec : records)
	{
		if (onlyPositiveStrand && !rec.id.strand()) continue;

		std::string contigSeq;
		for (size_t c = 0; c < rec.sequence.length(); c += FASTA_SLICE)
		{
			contigSeq += rec.sequence.substr(c, FASTA_SLICE).str() + "\n";
		}
		std::string header = onlyPositiveStrand ? 
							 ">" + rec.description.substr(1) + "\n":
							 ">" + rec.description + "\n";
		fwrite(header.data(), sizeof(header.data()[0]), 
			   header.size(), fout);
		fwrite(contigSeq.data(), sizeof(contigSeq.data()[0]), 
			   contigSeq.size(), fout);
	}
}

// By loadFrom File, each line was read and 
// (std::unordered_map<std::string, FastaRecord::Id> _nameIndex and vector<FastaRecord> _seqIndex) were filled
void SequenceContainer::buildPositionIndex()				// By this, _sequenceOffsets and _offsetsHint were filled
{
	Logger::get().debug() << "Building positional index";
	size_t offset = 0;
	_sequenceOffsets.reserve(_seqIndex.size());
	for (const auto& seq : _seqIndex)
	{
		_sequenceOffsets.push_back({offset, seq.sequence.length()});
		offset += seq.sequence.length();
	}
	_sequenceOffsets.push_back({offset, 0});
	if (offset == 0) return;

	_offsetsHint.reserve(offset / CHUNK + 1);
	size_t idx = 0;
	for (size_t i = 0; i <= (offset - 1) / CHUNK; ++i)
	{
		while (i * CHUNK >= _sequenceOffsets[idx + 1].offset) ++idx;
		//size_t newIdx = std::upper_bound(_sequenceOffsets.begin(), 
		//							     _sequenceOffsets.end(),
		//							     i * CHUNK) - _sequenceOffsets.begin();
		//Logger::get().debug() << idx << " " << newIdx;
		//assert(idx == newIdx);
		_offsetsHint.push_back(idx);
	}

	Logger::get().debug() << "Total sequence: " << offset / 2 << " bp";
	if (offset >= MAX_SEQUENCE)
	{
		Logger::get().error() << "Maximum sequence limit reached ("
			<< MAX_SEQUENCE / 2 << ")";
		throw std::runtime_error("Input overflow");
	}
}
