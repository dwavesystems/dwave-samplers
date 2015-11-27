#ifndef __ERRORS_H__
#define __ERRORS_H__

#include <exception>
#include <string>

class Errors: public std::exception{
public:
    std::string msg;
    Errors(std::string message): msg(message){};
    virtual const char* what() const throw(){
        return msg.c_str();
    }
    virtual ~Errors() throw(){};
};

#endif
