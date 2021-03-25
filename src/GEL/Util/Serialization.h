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

namespace Util {
    class Serialization
    {
        std::fstream fs;
        
    public:
        Serialization(const std::string& file_name, std::ios_base::openmode mode) {
            fs.open(file_name, std::ios_base::binary | mode);
        }
        
        template<typename T>
        void write(const T& blob) {
            fs.write(reinterpret_cast<const char*>(&blob), sizeof(T));
        }

        template<typename T>
        void read(T& blob) {
            fs.read(reinterpret_cast<char*>(&blob), sizeof(T));
        }
        
        template<typename T>
        void write(const std::vector<T>& vec) {
            size_t len = vec.size();
            const char* ptr = reinterpret_cast<const char*>(&vec[0]);
            fs.write(reinterpret_cast<const char*>(&len), sizeof(size_t));
            fs.write(ptr, len*sizeof(T));
        }
        
        template<typename T>
        void read(std::vector<T>& vec) {
            size_t len;
            fs.read(reinterpret_cast<char*>(&len), sizeof(size_t));
            vec.resize(len);
            fs.read(reinterpret_cast<char*>(&vec[0]), len*sizeof(T));
        }
        
        void write(const std::vector<bool>& vec) {
            size_t len = vec.size();
            std::vector<unsigned char> ucv(len);
            for(size_t i=0; i<len; ++i)
                ucv[i] = vec[i];
            write(ucv);
        }
        
        void read(std::vector<bool>& vec) {
            size_t len;
            fs.read(reinterpret_cast<char*>(&len), sizeof(size_t));
            std::vector<unsigned char> ucv(len);
            fs.read(reinterpret_cast<char*>(&ucv[0]), len*sizeof(unsigned char));
            vec.resize(len);
            for(size_t i=0; i<len; ++i)
                vec[i] = ucv[i];

        }

        
    };
}
#endif /* Serialization_hpp */
