#ifndef ERROR
#define ERROR

#include <string>
#include <stdio.h>

//namespace std{
class Error{

void throw_severe(std::string, int, bool);
void throw_minor(std::string, int, bool);

public:
	Error (std::string, int, bool); //general error contructor
	Error (std::string, int); //external-only error constructor
	Error (); //null constructor
	Error (const Error &); //copy
	Error operator = (const Error &);	//assignment


};
//}

#endif
