cmake_minimum_required(VERSION 3.10)
project(phasesplit C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(phasesplit
    src/main.c
    src/wav.c
    src/fft.c
    src/split.c
    src/dsp.c
    src/cpu.c
    src/thread.c
    src/fx.c
    src/upmix.c
    src/selftest.c
)
target_include_directories(phasesplit PRIVATE src)

# No -march or /arch switch anywhere on purpose. The wide paths are picked at
# run time after asking the processor what it has, so one binary starts on
# anything from an old SSE2 machine upwards. Building with -march=native would
# produce a faster binary that crashes on someone else's computer.
if(MSVC)
    target_compile_options(phasesplit PRIVATE /W4 /WX /O2 /fp:fast)
else()
    target_compile_options(phasesplit PRIVATE -Wall -Wextra -Werror -O2)
    # Reassociating floating point is fine for the spectral work and not fine
    # for the filters, which feed their own output back in. fx.c asks for
    # strict evaluation itself on MSVC; here the flag has to be given.
    set_source_files_properties(src/fx.c PROPERTIES COMPILE_OPTIONS "-ffp-contract=off")
    find_library(M_LIB m)
    if(M_LIB)
        target_link_libraries(phasesplit PRIVATE ${M_LIB})
    endif()
endif()

find_package(Threads REQUIRED)
target_link_libraries(phasesplit PRIVATE Threads::Threads)

install(TARGETS phasesplit RUNTIME DESTINATION bin)

enable_testing()
add_test(NAME selftest COMMAND phasesplit --selftest)
