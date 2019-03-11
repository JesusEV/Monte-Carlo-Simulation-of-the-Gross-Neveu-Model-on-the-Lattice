#ifndef COND_OS_STR_H
#define COND_OS_STR_H

#include <ostream>

class ConditionalOStream
{
    private:

        std::ostream &output_stream;
        bool active_flag = false;

    public:

        ConditionalOStream(std::ostream &stream, const bool active)
            : output_stream{stream}, active_flag{active}
        {}

        bool is_active() const noexcept { return active_flag; }

        void set_condition(bool flag) noexcept { active_flag = flag; }

        inline std::ostream& get_stream() const { return output_stream; }

        template <class T>
        inline const ConditionalOStream& operator<<(const T &t) const
        {
            if (active_flag == true) {output_stream << t;}
            return *this;
        }

        inline const ConditionalOStream& operator<<(std::ostream &(*p)(std::ostream &)) const
        {
            if (active_flag == true) {output_stream << p;}
            return *this;
        }
}; 
 
 #endif