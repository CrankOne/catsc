# Generic C/C++ implementation of cellular automata for track finding

This code is based on the algorithm of cellular automaton evolution
proposed for track finding by I.Abt et al [1].

For usage examples refer to `main-gpl-*.cc` files in the same directory.

## Design

Core functions are written in pure C to gain some performance wherever it
matters. Yet, note that most contributing factor (up to 70%) is still the
geometrical filter.

A C++ interface is provided in a form of single class with geometrical filter
and track candidates collector being an interfacing classes.

## Build

Standard CMake build procedure for ou-of-source build:

    $ cd build
    $ cmake ..
    $ make -j4 install

Note, that user might be interested in customizing type for hit identification
and/or floating point type to describe spatial coordinates -- for that purpose
once can set CMake variables `CATS_HIT_ID_TYPE` and `CATS_COORDINATE_TYPE`
with `-D...=...` respictively.

## Usage in user code

For CMake-based projects:

    find_package(catsc)
    # ...
    target_link_libraries(${na64std_LIB} PUBLIC catsc::catsc )

For other types of build system a `pkg-config` file is exported by given user
prefix.

# References

[1] Abt, I.; Emeliyanov, D.; Kisel, I.; Masciocchi, S.; CATS: a cellular
automaton for tracking in silicon for the HERA-B vertex detector // Nuclear
Instruments and Methods in Physics Research A 489 (2002) 389–405

