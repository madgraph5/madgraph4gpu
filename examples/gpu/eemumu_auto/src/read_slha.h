#ifndef READ_SLHA_H
#define READ_SLHA_H

#include <map>
#include <string>
#include <sstream>
#include <vector>
#include "mg5Complex.h"

class SLHABlock
{
  public:
    SLHABlock(std::string name = ""){_name = name;}
    ~SLHABlock(){}

    void set_entry(std::vector<int> indices, floa_t value);
    floa_t get_entry(std::vector<int> indices, floa_t def_val = 0);
    void set_name(std::string name) {_name = name;}
    std::string get_name(){return _name;}
    int get_indices() { return _indices;}

  private:
    std::string _name;
    std::map<std::vector<int>, double> _entries;
    unsigned int _indices;
};

class SLHAReader
{
  public:
    SLHAReader(std::string file_name = "")
	{if(file_name != "") read_slha_file(file_name);}

    void read_slha_file(std::string file_name);
    floa_t get_block_entry(std::string block_name, std::vector<int> indices,
			   floa_t def_val = 0);
    floa_t get_block_entry(std::string block_name, int index,
			   floa_t def_val = 0);
    void set_block_entry(std::string block_name, std::vector<int> indices,
			   floa_t value);
    void set_block_entry(std::string block_name, int index,
			   floa_t value);
  private:
    std::map<std::string, SLHABlock> _blocks;
};

#endif
