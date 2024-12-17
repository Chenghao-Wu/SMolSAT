/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/


#ifndef TOKENIZE
#define TOKENIZE
#include <string>
#include <vector>
#include "vector_map.h"

#define ARGMAX 10000

//namespace std{

class Tokenize
{
    std::vector <std::string> tokens;
    std::Vector_Map <std::string,std::string> flags;
    bool not_empty;
    void parse_token(std::string line, size_t token_start);
    void parse_token(std::string line, size_t token_start, size_t token_end);
    std::string flagmarker;
    bool flagcheck;

public:
    Tokenize(){flagcheck=false;flagmarker="";};			//default constructor just sets flag checking to off
    Tokenize(std::string line);		//constructor to tokenize initial line
    Tokenize(const Tokenize & copy);	//copy constructor
    Tokenize operator=(const Tokenize & copy);	//equality operator
    //void operator()(std::string line);	//tokenize new line
    std::vector <std::string> operator()(std::string line);
    
    void setflagmarker(std::string marker){flagmarker=marker;flagcheck=true;};	//set the flagmarker and set to check for flags
    void setflagmarker(){flagmarker="";flagcheck=false;};	//clear the flagmarker and set to not check for flags
    
    std::string operator()(int tokenindex){return tokens[tokenindex];};	//return token corresponding to tokenindex
    std::string operator[](std::string flagid);
    int operator()(std::string* tokenarray,int maxsize);	//returns number of tokens and sets tokenarray* to be an array of tokens up to size maxsize
    int operator()(std::string line, std::string* tokenarray, int maxsize=ARGMAX);	//tokenizes new line, returns number of tokens, and sets tokenarray to contain tokens
	
    
    //int operator(std::string line,std::string* tokenlist);
    int count(){return tokens.size();};		//returns number of tokens
    bool isflagged(std::string flagid){return flags.count(flagid);};		//returns a 1 if flag given by flag id is in line; 0 otherwise.
    
    
    bool in_string_array(std::string target);	//check if a token is in the line
    int find_in_string_array(std::string target);	//return the position of a token if in the line
};

//}
  
#endif


