%module odsl
%{
#include "interface/common_interface.hpp"
#include "interface/omap_interface.hpp"
#include "interface/par_omap_interface.hpp"
#include "interface/recoram_interface.hpp"
%}

/* Let's just grab the original header file here */
%include "interface/common_interface.hpp"
%include "interface/omap_interface.hpp"
%include "interface/par_omap_interface.hpp"
%include "interface/recoram_interface.hpp"