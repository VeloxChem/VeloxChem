//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef FmtType_hpp
#define FmtType_hpp

/**
 Enumerate class fmt:

 Defines supported formatting keys:
 fmt::end    - the end of output line
 fmt::left   - the left alignment of output line
 fmt::center - the center alignment of output line
 fmt::right  - the right alignment of output line
 fmt::blank  - the empty output line terminated by fmt::end
 fmt::title  - the title style of output line
 fmt::tsep   - the separator output line terminated by fmt::end
 fmt::cerror - the critial error style of output line
 fmt::into   - the information message style of output line
 fmt::error  - the syntax error style of output line
 fmt::header - the header style of output line
 */
enum class fmt
{
    end,
    left,
    center,
    right,
    blank,
    title,
    tsep,
    cerror,
    info,
    error,
    header
};

#endif /* FmtType_hpp */
