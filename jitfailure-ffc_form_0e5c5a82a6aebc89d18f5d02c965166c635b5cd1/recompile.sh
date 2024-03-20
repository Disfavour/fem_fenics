#!/bin/bash
# Execute this file to recompile locally
c++ -Wall -shared -fPIC -std=c++11 -O2 -I/usr/lib/python3/dist-packages/ffc/backends/ufc -I/home/disfavour/.cache/dijitso/include ffc_form_0e5c5a82a6aebc89d18f5d02c965166c635b5cd1.cpp -L/home/disfavour/.cache/dijitso/lib -Wl,-rpath,/home/disfavour/.cache/dijitso/lib -ldijitso-ffc_element_0aceea476c4466c38bcd4b3da92b7c720101a8ac -ldijitso-ffc_element_f15c62f5d90fd349915de9977c93d95ae6a6e4ca -ldijitso-ffc_coordinate_mapping_7b2d1da84570d09b9efefe42fa819358fb99594c -olibdijitso-ffc_form_0e5c5a82a6aebc89d18f5d02c965166c635b5cd1.so