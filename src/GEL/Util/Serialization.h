//
//  Serialization.hpp
//  GEL
//
//  Created by Jakob Andreas Bærentzen on 24/01/2018.
//  Copyright © 2018 J. Andreas Bærentzen. All rights reserved.
//

#ifndef Serialization_hpp
#define Serialization_hpp

#include <vector>
#include <string>
#include <fstream>

class Serialization
{
    std::fstream fs;
    
public:
    Serialization(const std::string& file_name, std::ios_base::openmode mode) {
        fs.open(file_name, std::ios_base::binary | mode);
    }
    
    void write(const char* blob, size_t len) {
        fs.write(blob, len);
    }
    void read(char* blob, size_t len) {
        fs.read(blob, len);
    }
    
    template<typename T>
    void write(const std::vector<T>& vec) {
        size_t len = vec.size();
        char* ptr = reinterpret_cast<char*>(&vec[0]);
        write(reinterpret_cast<char*>(&len), sizeof(size_t));
        write(ptr,len * sizeof(T));
    }

    template<typename T>
    void read(std::vector<T>& vec) {
        size_t len;
        read(reinterpret_cast<char*>(&len), sizeof(size_t));
        vec.resize(len);
        read(reinterpret_cast<char*>(&vec[0]), len*sizeot(T));
    }

};

#endif /* Serialization_hpp */
