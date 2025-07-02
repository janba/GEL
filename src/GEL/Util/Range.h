//
//  NumericRange.h
//  GEL
//
//  Created by J. Andreas Bærentzen on 17/02/2016.
//  Copyright © 2016 J. Andreas Bærentzen. All rights reserved.
//

#ifndef NumericRange_h
#define NumericRange_h

#include <ranges>

namespace Util {
// FIXME: this just duplicates functionality of std::ranges::iota_view, consider deprecating it
    class Range {
    public:
        class iterator {
            friend class Range;
        public:
            // typedefs to accommodate stl compliance
            typedef ptrdiff_t difference_type;
            typedef std::forward_iterator_tag iterator_category;
            typedef long int value_type;
            typedef value_type reference;
            typedef value_type* pointer;

            value_type operator *() const { return i_; }
            iterator &operator ++() { ++i_; return *this; }
            iterator operator ++(int) { iterator copy(*this); ++i_; return copy; }

            bool operator ==(const iterator &other) const { return i_ == other.i_; }
            bool operator !=(const iterator &other) const { return i_ != other.i_; }

            // required for std::viewable_range
            iterator() = default;
        protected:
            iterator(long int start) : i_ (start) { }

        private:
            unsigned long i_;
        };

        iterator begin() const { return begin_; }
        iterator end() const { return end_; }
        Range(long int  begin, long int end) : begin_(begin), end_(end) {}
    private:
        iterator begin_;
        iterator end_;
    };
    static_assert(std::input_or_output_iterator<Range::iterator>);
    static_assert(std::ranges::viewable_range<const Range&>);
}

#endif /* NumericRange_h */
