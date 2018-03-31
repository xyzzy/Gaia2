cmake_minimum_required(VERSION 3.9)
project(gaia)

set(CMAKE_CXX_STANDARD 98)
add_executable(gaia2 gaia2.cc)
add_executable(merge merge.cc)
add_executable(simulTame4 simulTame4.cc)

target_include_directories(gaia2 PRIVATE /pub/disk42/untangle/jansson-2.10/destdir/include)

# compiled with:
# g++ simulTame4.cc -o simulTame4 /usr/lib/libgd.a /usr/lib/libpng.a /usr/lib/libjpeg.a  /usr/lib/libz.a -lpthread -O3
# I use many experinental libraries. I link to statics because I know they are safe and correct.
# run:
# ./simulTame4 --seed=00001 images/heen-0000d.%02d.png weer-0000.png heen-0000.png 25
# it can take some 10 minutes per frame to calculate