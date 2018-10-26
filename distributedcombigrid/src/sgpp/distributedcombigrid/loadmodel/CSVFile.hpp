
#ifndef CSVFILE_HPP_ 
#define CSVFILE_HPP_ 

#ifndef USE_HDF5
// csv file class, adapted from
// https://gist.github.com/rudolfovich/f250900f1a833e715260a66c87369d15

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include <boost/tokenizer.hpp>

class csvfile;

inline static csvfile& endrow(csvfile& file);
inline static csvfile& flush(csvfile& file);

class csvfile
{
    std::string filename_;
    std::ofstream fs_;
    bool is_first_;
    const char separator_;
    const std::string escape_seq_;
    const std::string special_chars_;
public:
    csvfile(const std::string filename, const char separator = ',')
        : filename_(filename)
        , fs_()
        , is_first_(true)
        , separator_(separator)
        , escape_seq_("\"")
        , special_chars_("\"")
    {
        fs_.exceptions(std::ios::failbit | std::ios::badbit);
        fs_.open(filename, std::ofstream::app);
    }

    ~csvfile()
    {
        flush();
        fs_.close();
    }

    //only to be used in test functions!
    void remove_file(){
        remove(filename_.c_str());
    }

    void flush()
    {
        fs_.flush();
    }

    void endrow()
    {
        fs_ << std::endl;
        is_first_ = true;
    }

    csvfile& operator << ( csvfile& (* val)(csvfile&))
    {
        return val(*this);
    }

    csvfile& operator << (const char * val)
    {
        return write(escape(val));
    }

    csvfile& operator << (const std::string & val)
    {
        return write(escape(val));
    }

    template<typename T>
    csvfile& operator << (const T& val)
    {
        return write(val);
    }

    //returns a vector of the contents of a column (zero-based indexing), as long as they can be cast to T
    template<typename T>
    std::vector<T> readColumn (size_t column) const
    {
        std::ifstream in(filename_.c_str());
        assert (in.is_open());

        std::string line;
        std::vector<T> vec;

        while (getline(in,line)){
            boost::tokenizer<boost::escaped_list_separator<char> > tok(line, boost::escaped_list_separator<char>('\\', separator_, '\"'));
            boost::tokenizer<boost::escaped_list_separator<char> >::iterator i(tok.begin());
            for (size_t j=0 ; j < column; ++j){
                ++i;
            }
            try{
                T t = std::strtol ((*i).c_str(),nullptr,10);
                vec.push_back(t);
            }
            catch (const std::invalid_argument& ia) {
            }
        }
        return vec;
    }

private:
    template<typename T>
    csvfile& write (const T& val)
    {
        if (!is_first_)
        {
            fs_ << separator_;
        }
        else
        {
            is_first_ = false;
        }
        fs_ << val;
        return *this;
    }

    std::string escape(const std::string & val)
    {
        std::ostringstream result;
        result << '"';
        std::string::size_type to, from = 0u, len = val.length();
        while (from < len &&
                std::string::npos != (to = val.find_first_of(special_chars_, from)))
        {
            result << val.substr(from, to - from) << escape_seq_ << val[to];
            from = to + 1;
        }
        result << val.substr(from) << '"';
        return result.str();
    }
};


inline static csvfile& endrow(csvfile& file)
{
    file.endrow();
    return file;
}

inline static csvfile& flush(csvfile& file)
{
    file.flush();
    return file;
}

#endif //def USE_HDF5

#endif /* CSVFILE_HPP_ */